"""The LED flip-detection network: a dilated 1-D convolution stack with slot readout.

Input is the normalised I/Q record ``[1, 2, N]`` plus side inputs added to the first hidden activation: the LED timing
cues (time since the latest onset, LED-on mask), three per-sample channels from the LED-window HMM, three per-sample
channels from its smoothed-fit variant (HMM-SF), and one record-level burst-intensity feature. Eleven residual blocks
(dilations 1 ... 1024, ~41 ms receptive field at 100 kHz) feed a shared head that reads out slots: the whole record is
tiled into background groups at a fixed density, plus a denser tiling of the 20 ms after every LED onset. Every slot
predicts a presence probability and a time inside its support. 9848 learned values in all.
"""
import math

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F

ON_TIME_S = 50e-6
OFF_TIME_S = .300
PERIOD_S = ON_TIME_S + OFF_TIME_S
HIDDEN = 16
DILATIONS = (1, 2, 4, 8, 16, 32, 64)
CONTEXT_DILATIONS = (128, 256)
CONTEXT2_DILATIONS = (512, 1024)
SIDE_CHANNELS = 3
GEOMETRY = ('gather_i0', 'gather_i1', 'gather_weight', 'support_start', 'support_width', 'anchor_logits')
N_PARAMETERS = 9848
WEIGHTS_FORMAT = 'qpd-led-flip-network/1'


def phase_channels(n_samples, first_pulse_time_s, sample_rate_hz=100000., *, device=None, dtype=torch.float32):
    """[1, 2, N] LED timing cues computed from the schedule alone.

    Channel 0 is the time since the most recent LED rise divided by the period: in [0, 1) after the first rise, and
    the (negative) time to the first rise divided by the period before it. Channel 1 marks [rise, rise + 50 us).
    """
    if isinstance(n_samples, bool) or int(n_samples) != n_samples or n_samples < 1:
        raise ValueError("n_samples must be a positive integer")
    first, rate = float(first_pulse_time_s), float(sample_rate_hz)
    if not math.isfinite(first) or first < 0 or not math.isfinite(rate) or rate <= 0:
        raise ValueError("First onset must be finite/nonnegative and sample rate finite/positive")
    if dtype not in (torch.float32, torch.float64):
        raise ValueError("Timing cues require float32 or float64")
    # ON and OFF are added in sample units separately, so the default boundaries are exactly 5 and 30005 samples.
    on_samples = ON_TIME_S * rate
    period_samples = on_samples + OFF_TIME_S * rate
    first_samples = first * rate
    if not all(math.isfinite(v) for v in (first_samples, on_samples, period_samples)):
        raise ValueError("Schedule sample coordinates must remain finite")
    relative = torch.arange(int(n_samples), device=device, dtype=torch.float64) - first_samples
    after = relative >= 0
    elapsed = torch.where(after, torch.remainder(relative, period_samples), relative)
    led_on = after & (elapsed < on_samples)
    phase = (elapsed / period_samples).to(dtype)
    # Casting an almost-complete cycle to float32 must not turn it into 1.
    upper = torch.nextafter(phase.new_tensor(1), phase.new_tensor(0))
    phase = torch.minimum(phase, upper)
    return torch.stack((phase, led_on.to(dtype)))[None]


def _block(dilation):
    return nn.Conv1d(HIDDEN, HIDDEN, 3, dilation=dilation, padding=dilation)


class LEDFlipNetwork(nn.Module):
    """The network. Build it with :func:`load_network`; set each record with :func:`configure` and the
    ``configure_*`` methods, then call :func:`forward`."""

    def __init__(self):
        super().__init__()
        self.hidden_channels = HIDDEN
        self.input_projection = nn.Conv1d(2, HIDDEN, 1)
        self.residual_blocks = nn.ModuleList([_block(d) for d in DILATIONS])
        self.shared_slot_head = nn.Linear(8 * HIDDEN, 8)
        self.cue_projection = nn.Conv1d(2, HIDDEN, 1, bias=False)
        self.hmm_projection = nn.Conv1d(SIDE_CHANNELS, HIDDEN, 1, bias=False)
        self.context_blocks = nn.ModuleList([_block(d) for d in CONTEXT_DILATIONS])
        self.context2_blocks = nn.ModuleList([_block(d) for d in CONTEXT2_DILATIONS])
        self.sf_projection = nn.Conv1d(SIDE_CHANNELS, HIDDEN, 1, bias=False)
        self.mu_projection = nn.Linear(1, HIDDEN, bias=False)
        for name in GEOMETRY:
            self.register_buffer(name, torch.empty(0), persistent=False)
        self.register_buffer('_phase_cues', torch.empty((1, 2, 0)), persistent=False)
        self.register_buffer('_hmm_channels', torch.empty((1, SIDE_CHANNELS, 0)), persistent=False)
        self.register_buffer('_sf_channels', torch.empty((1, SIDE_CHANNELS, 0)), persistent=False)
        self._mu_value = None
        self._geometry_5s = None

    def configure_phase(self, n_samples, first_pulse_time_s, sample_rate_hz=100000.):
        weight = self.input_projection.weight
        self._phase_cues = phase_channels(n_samples, first_pulse_time_s, sample_rate_hz, device=weight.device, dtype=weight.dtype)
        return self

    def configure_hmm(self, channels):
        """[1, 3, N] LED-window HMM channels (swing, burst, confidence)."""
        weight = self.input_projection.weight
        if channels.ndim != 3 or channels.shape[:2] != (1, SIDE_CHANNELS):
            raise ValueError('Expected [1,3,N] HMM channels')
        self._hmm_channels = channels.to(device=weight.device, dtype=weight.dtype)
        return self

    def configure_sf(self, channels):
        """[1, 3, N] HMM-SF channels (contrast, gate, posterior swing)."""
        weight = self.input_projection.weight
        if channels.ndim != 3 or channels.shape[:2] != (1, SIDE_CHANNELS):
            raise ValueError('Expected [1,3,N] HMM-SF channels')
        self._sf_channels = channels.to(device=weight.device, dtype=weight.dtype)
        return self

    def configure_mu(self, value):
        """value: None (feature off) or the record's log1p(mu_hat)/log1p(10)."""
        if value is None:
            self._mu_value = None
            return self
        value = float(value)
        if not value == value or abs(value) == float('inf'):
            raise ValueError('mu input must be finite')
        weight = self.input_projection.weight
        self._mu_value = torch.tensor(value, device=weight.device, dtype=weight.dtype)
        return self

    def forward(self, x):
        n = x.shape[-1]
        if (x.ndim != 3 or x.shape[1] != 2 or self._phase_cues.shape[-1] != n or self._hmm_channels.shape[-1] != n
                or self._sf_channels.shape[-1] != n):
            raise ValueError('Expected [B,2,N]; configure_phase/configure_hmm/configure_sf must match this record')
        h = (self.input_projection(x.float())+self.cue_projection(self._phase_cues)+self.hmm_projection(self._hmm_channels)
             + self.sf_projection(self._sf_channels))
        if self._mu_value is not None:
            h = h+(self.mu_projection.weight[:, 0]*self._mu_value)[None, :, None]
        for block in list(self.residual_blocks)+list(self.context_blocks)+list(self.context2_blocks):
            h = h+F.silu(block(h))
        left = h[:, :, self.gather_i0]
        right = h[:, :, self.gather_i1]
        weight = self.gather_weight.unsqueeze(0).unsqueeze(0)
        gathered = (1.0-weight)*left+weight*right
        group_features = gathered.permute(0, 2, 3, 1).reshape(gathered.shape[0], gathered.shape[2], 8*self.hidden_channels)
        head_output = self.shared_slot_head(group_features)
        presence = torch.sigmoid(head_output[..., :4])
        fraction = torch.sigmoid(self.anchor_logits.reshape(1, -1, 4)+head_output[..., 4:])
        sample_time = self.support_start.reshape(1, -1, 4)+self.support_width.reshape(1, -1, 4)*fraction
        normalized_time = sample_time/500000.0
        return torch.stack([normalized_time.reshape(normalized_time.shape[0], -1),
                            presence.reshape(presence.shape[0], -1)], dim=-1).float()


def load_network(path):
    """Build the network and load its learned values (plus the fixed 5 s slot geometry) from ``path``."""
    obj = torch.load(path, map_location='cpu', weights_only=True)
    if not isinstance(obj, dict) or obj.get('format') != WEIGHTS_FORMAT:
        raise ValueError(f'{path}: not a qpd LED flip-network weights file')
    model = LEDFlipNetwork()
    params = dict(model.named_parameters())
    state = obj['parameters']
    if set(state) != set(params) or sum(v.numel() for v in state.values()) != N_PARAMETERS:
        raise ValueError(f'{path}: expected exactly the {N_PARAMETERS} learned values of the network')
    with torch.no_grad():
        for k, p in params.items():
            if p.shape != state[k].shape:
                raise ValueError(f'{path}: {k} has shape {tuple(state[k].shape)}, expected {tuple(p.shape)}')
            p.copy_(state[k])
    model._geometry_5s = {k: (v.to(torch.int64) if k in ('gather_i0', 'gather_i1') else v.clone())
                          for k, v in obj['geometry_5s'].items()}
    for k, v in model._geometry_5s.items():
        setattr(model, k, v.clone())
    return model.eval()


def configure(model, n_samples, first_pulse_time_s, sample_rate_hz=100000.):
    """Set the LED cues and tile the record into slots.

    The whole record is tiled at the trained density, and each available 20 ms post-onset tail at the denser LED
    density. A 5 s record whose first pulse is at exactly 50 ms uses the fixed geometry stored with the weights.
    """
    model.configure_phase(n_samples, first_pulse_time_s, sample_rate_hz)
    dev = model.support_start.device; dtype = model.support_start.dtype
    n = int(n_samples); end = float(np.nextafter(np.float32(n), np.float32(0)))
    if n == 500000 and abs(first_pulse_time_s-.05) < 1e-12:
        for name, value in model._geometry_5s.items():
            setattr(model, name, value.clone())
        return
    groups = max(1, round(4096*n/500000))
    start = torch.arange(groups, device=dev, dtype=dtype)*(end/groups)
    width = torch.full_like(start, end/groups)
    rises = np.arange(float(first_pulse_time_s), n/sample_rate_hz, .30005)
    starts = [start]; widths = [width]
    for rise in rises:
        lo = rise*sample_rate_hz; hi = min(lo+2000., end)
        if hi <= lo:
            continue
        count = max(1, round((hi-lo)*241/2000.))
        starts.append(lo+torch.arange(count, device=dev, dtype=dtype)*(hi-lo)/count)
        widths.append(torch.full((count,), (hi-lo)/count, device=dev, dtype=dtype))
    start, width = torch.cat(starts), torch.cat(widths)
    q = (torch.arange(8, device=dev, dtype=dtype)+.5)/8
    coordinates = (start[:, None]+width[:, None]*q).clamp(0, n-1)
    left = coordinates.floor().long()
    values = dict(gather_i0=left, gather_i1=(left+1).clamp_max(n-1), gather_weight=coordinates-left,
                  support_start=start.repeat_interleave(4), support_width=width.repeat_interleave(4),
                  anchor_logits=torch.logit(torch.tensor([.125, .375, .625, .875], device=dev, dtype=dtype)).repeat(len(start)))
    for name, value in values.items():
        setattr(model, name, value)


def forward(model, x):
    """Run the network; returns [B, slots, 2] = (time as a fraction of the record, presence probability)."""
    out = model(x)
    return torch.stack((out[..., 0]*(500000/x.shape[-1]), out[..., 1]), -1)

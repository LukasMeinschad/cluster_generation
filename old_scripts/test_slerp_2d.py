# -*- coding: utf-8 -*-
"""
Simple 2D SLERP visualization.
All rotations are around the z-axis, so we can plot everything in the xy-plane.

Test 1: Rotate a vector from 0 to 120 deg - SLERP vs LERP on the unit circle
Test 2: Angle vs t - shows SLERP gives constant angular velocity, LERP does not
Test 3: Arrow fan - SLERP interpolated arrows at equal t-steps are evenly spaced
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "modules"))

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch
from transformations import Quaternion

# ============ Setup ============
# Two rotations around z-axis: 0° and 120°
angle_start = 0.0
angle_end = 120.0

q_start = Quaternion.from_axis_angle(np.array([0, 0, 1]), np.radians(angle_start))
q_end   = Quaternion.from_axis_angle(np.array([0, 0, 1]), np.radians(angle_end))

N = 50
ts = np.linspace(0, 1, N)
test_vec = np.array([1.0, 0.0, 0.0])  # unit vector along x


def lerp_quat(q1, q2, t):
    """Naive linear interpolation + normalize (for comparison)."""
    return Quaternion(
        q1.w + t * (q2.w - q1.w),
        q1.x + t * (q2.x - q1.x),
        q1.y + t * (q2.y - q1.y),
        q1.z + t * (q2.z - q1.z),
    ).normalize()


# ============ Compute paths ============
slerp_xy = []
lerp_xy  = []
slerp_angles = []
lerp_angles  = []

for t in ts:
    qs = Quaternion.slerp(q_start, q_end, t)
    ql = lerp_quat(q_start, q_end, t)

    ps = qs.rotate_point(test_vec)
    pl = ql.rotate_point(test_vec)

    slerp_xy.append(ps[:2])
    lerp_xy.append(pl[:2])
    slerp_angles.append(np.degrees(np.arctan2(ps[1], ps[0])))
    lerp_angles.append(np.degrees(np.arctan2(pl[1], pl[0])))

slerp_xy = np.array(slerp_xy)
lerp_xy  = np.array(lerp_xy)

# ============ Figure ============
fig, axes = plt.subplots(1, 3, figsize=(18, 5.5))
fig.suptitle("SLERP vs LERP -- 2D Rotation Test (0 deg -> 120 deg around z)", fontsize=14, fontweight='bold')

# ---------- Plot 1: Paths on unit circle ----------
ax = axes[0]
theta_circle = np.linspace(0, 2 * np.pi, 300)
ax.plot(np.cos(theta_circle), np.sin(theta_circle), 'k-', lw=0.4, alpha=0.3)

ax.plot(slerp_xy[:, 0], slerp_xy[:, 1], 'b-o', ms=3, lw=2, label='SLERP', zorder=3)
ax.plot(lerp_xy[:, 0],  lerp_xy[:, 1],  'r--s', ms=2.5, lw=1.5, label='LERP (renorm.)', zorder=2)

# Endpoints
ax.scatter(*slerp_xy[0],  color='green', s=120, zorder=5, edgecolors='k', label='Start (0 deg)')
ax.scatter(*slerp_xy[-1], color='orange', s=120, zorder=5, edgecolors='k', label='End (120 deg)')

# t-labels on SLERP path
for i in range(0, N, N // 5):
    ax.annotate(f't={ts[i]:.1f}', xy=slerp_xy[i], fontsize=7, color='blue',
                xytext=(5, 5), textcoords='offset points')

ax.set_aspect('equal')
ax.set_xlim(-1.3, 1.3)
ax.set_ylim(-0.4, 1.4)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('Paths on Unit Circle')
ax.legend(fontsize=8, loc='lower left')
ax.grid(True, alpha=0.2)

# ---------- Plot 2: Angle vs t ----------
ax = axes[1]
ax.plot(ts, slerp_angles, 'b-', lw=2.5, label='SLERP')
ax.plot(ts, lerp_angles,  'r--', lw=2, label='LERP (renorm.)')
ax.plot(ts, ts * 120, 'k:', lw=1.5, label='Ideal linear (120 deg * t)')

ax.set_xlabel('t', fontsize=12)
ax.set_ylabel('Angle (degrees)', fontsize=12)
ax.set_title('Angular Velocity Check')
ax.legend(fontsize=9)
ax.grid(True, alpha=0.3)

# Residual inset
inset = ax.inset_axes([0.55, 0.1, 0.4, 0.35])
inset.plot(ts, np.array(slerp_angles) - ts * 120, 'b-', lw=1.5, label='SLERP error')
inset.plot(ts, np.array(lerp_angles) - ts * 120, 'r--', lw=1.5, label='LERP error')
inset.set_ylabel('delta angle (deg)', fontsize=7)
inset.set_xlabel('t', fontsize=7)
inset.tick_params(labelsize=6)
inset.legend(fontsize=6)
inset.grid(True, alpha=0.2)
inset.set_title('Deviation from ideal', fontsize=7)

# ---------- Plot 3: Arrow fan ----------
ax = axes[2]
ax.plot(np.cos(theta_circle), np.sin(theta_circle), 'k-', lw=0.4, alpha=0.3)

n_arrows = 12
arrow_ts = np.linspace(0, 1, n_arrows)
colors = plt.cm.viridis(np.linspace(0, 1, n_arrows))

for i, t in enumerate(arrow_ts):
    qs = Quaternion.slerp(q_start, q_end, t)
    p = qs.rotate_point(test_vec)
    ax.annotate('', xy=(p[0], p[1]), xytext=(0, 0),
                arrowprops=dict(arrowstyle='->', color=colors[i], lw=2))
    # Small text label
    ax.text(p[0] * 1.12, p[1] * 1.12, f'{t:.2f}', fontsize=7,
            ha='center', va='center', color=colors[i])

ax.set_aspect('equal')
ax.set_xlim(-1.4, 1.4)
ax.set_ylim(-0.5, 1.5)
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_title('SLERP Arrow Fan (equal t-spacing)')
ax.grid(True, alpha=0.2)

# Colorbar
sm = plt.cm.ScalarMappable(cmap='viridis', norm=plt.Normalize(0, 1))
sm.set_array([])
cbar = plt.colorbar(sm, ax=ax, shrink=0.7, label='t parameter')

plt.tight_layout()
save_path = os.path.join(os.path.dirname(__file__), "modules", "figures", "slerp_2d_test.png")
plt.savefig(save_path, dpi=150, bbox_inches='tight')
print(f"Saved to {save_path}")
plt.show()

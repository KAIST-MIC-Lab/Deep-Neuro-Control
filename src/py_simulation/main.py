
import numpy as np
from Env import Env
import matplotlib.pyplot as plt

# ============================================
#         SIMULATION VARIABLES
# ============================================
RESULT_SAVE = False
FIGURE_PLOT = True
FIGURE_SAVE = False

T = 10.0            # simulation time
ctrl_dt = 1e-4      # control time step
dt = ctrl_dt * 1e-1 # time step
rpt_dt = 0.1        # report time step (on console)
t = np.arange(0, T, dt)

# ============================================
#         LOCAL FUNCTIONS
# ============================================
def ref_gen(t, x):
    r = np.array([
        np.sin(t),
        np.cos(t)
    ]).reshape(-1, 1)

    return r

# ============================================
#         MAIN FUNCTION
# ============================================
if __name__ == "__main__":
    env = Env(dt)
    x = np.array([[10], [0]])
    env.reset(x)

    x_hist = []
    u_hist = []
    r_hist = []

    K = np.array([
        10, 15
    ]).reshape(1, -1)
    
    print("Simulation Start")

    for t_idx in range(len(t)):

        if (t_idx % (ctrl_dt / dt)) == 0:
            x = env.x
            xd = ref_gen(t[t_idx], env.x)
            e = x - xd
            u = -np.dot(K, e)

        env.step(u)

        x_hist.append(env.x)
        u_hist.append(u)
        r_hist.append(xd)

        if (t_idx % (rpt_dt / dt)) == 0:
            print(f"Simulation Time: {t_idx * dt:.3f} sec")

    print("Simulation Done")

    if FIGURE_PLOT:
        font_size = 16;
        font_size_ticl = 12;
        lgd_size = 12;
        line_width = 2;
        fig_height = 5; 
        fig_width = 6;

        fontdict={'fontname': 'Times New Roman',
            'fontsize': font_size,
            'fontweight': 'bold'}  # 'heavy', 'light', 'ultrabold', 'ultralight'

        plt.rcParams['text.usetex'] = True
        plt.rcParams['font.family'] = 'Times New Roman'
        plt.rcParams['font.size'] = font_size_ticl 

        plt.figure(1, figsize=(fig_width, fig_height))
        plt.plot(t, [x[0] for x in x_hist], label="$x_1$", color="blue", linewidth=line_width)
        plt.plot(t, [xd[0] for xd in r_hist], label="$r_1$", color="red", linewidth=line_width)
        plt.title("Tracking Result of $x_1$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$x_1$',  fontdict=fontdict);

        plt.figure(2, figsize=(fig_width, fig_height))
        plt.plot(t, [x[1] for x in x_hist], label="$x_2$", color="blue", linewidth=line_width)
        plt.plot(t, [xd[1] for xd in r_hist], label="$r_2$", color="red", linewidth=line_width)
        plt.title("Tracking Result of $x_2$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$x_2$',  fontdict=fontdict);

        plt.figure(3, figsize=(fig_width, fig_height))
        plt.plot(t, [u[0] for u in u_hist], label="$u$", color="blue", linewidth=line_width)
        plt.title("Control Input", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$u$',  fontdict=fontdict);
        
        plt.show()

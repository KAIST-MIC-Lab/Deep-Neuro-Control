
import numpy as np
from EnvSimple import Env
import matplotlib.pyplot as plt
from LMI_controller import CVSTEM
from NCM import NCM
import pickle

# ============================================
#         SIMULATION VARIABLES
# ============================================
RESULT_SAVE = 0
FIGURE_PLOT = 1
FIGURE_SAVE = 0

NN_PATH = "train_src/model.pth"

T = 4.0            # simulation time
ctrl_dt = 1e-2      # control time step
dt = ctrl_dt * 1e0 # time step
rpt_dt = 1        # report time step (on console)
t = np.arange(0, T, dt)

# ============================================
#         LOCAL FUNCTIONS
# ============================================
def ref_input(t, x):
    ud = np.array([
        np.sin(t)+1,
        -np.cos(t)+1
    ]).reshape(-1, 1) * 10

    return ud

def simulate(x0, xd0):
    env = Env(dt)
    env_ref = Env(dt)

    x = x0
    # x = np.array([
    #     0.0, 
    #     0.0, 
    # ], dtype="f").reshape(2,1);
    
    env.reset(x)

    xd0 = xd0
    # xd0 = np.array([
    #     10.0, 
    #     -10.0, 
    # ], dtype="f").reshape(2,1);

    env_ref.reset(xd0)

    x_hist = []
    u_hist = []
    r_hist = []
    ud_hist = []
    M_hist = []

    # ctrl_NCM = NCM(NN_PATH)
    ctrl_CVSTEM = CVSTEM(ctrl_dt)
    
    print("Simulation Start")

    for t_idx in range(len(t)):
        ud = ref_input(t[t_idx], env_ref.x)
        env_ref.step(ud)

        if (t_idx % (ctrl_dt / dt)) == 0:
            x = env.x
            xd = env_ref.x

            # u, M = ctrl_NCM.getControl(x, xd, ud)
            # u_CVSTEM, M_CVSTEM = ctrl_CVSTEM.getControl(x, xd, ud)
            u, M = ctrl_CVSTEM.getControl(x, xd, ud)

            # print(f"norm of M: {np.linalg.norm(M)}")
            # print(f"norm of M_CVSTEM: {np.linalg.norm(M_CVSTEM)}")

            # M = np.zeros((2, 2))
            # u = ud
            # u = np.zeros((2, 1), dtype="f")

        env.step(u)

        x_hist.append(env.x)
        u_hist.append(u)
        r_hist.append(xd)
        ud_hist.append(ud)
        M_hist.append(np.linalg.norm(M,2))

        if (t_idx % (rpt_dt / dt)) == 0:
            print(f"Simulation Time: {t_idx * dt:.3f} sec")

    print("Simulation Done")

    return x_hist, u_hist, r_hist, ud_hist, M_hist
    
def result_plot():
    # =================================
        #         PLOT CONFIGURATION
        # =================================
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

        # =================================
        #   FIG 1: Error (Alpha)
        # =================================
        plt.subplot(241)
        plt.plot(t, [x[0] - xd[0] for x, xd in zip(x_hist, r_hist)], label="$e_{\\alpha}$", color="blue", linewidth=line_width)
        plt.title("Tracking Error of $e_{\\alpha}$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$e_{\\alpha}$ / A',  fontdict=fontdict);

        # =================================
        #   FIG 2: Error (Beta)
        # =================================
        plt.subplot(242)
        plt.plot(t, [x[1] - xd[1] for x, xd in zip(x_hist, r_hist)], label="$e_{\\beta}$", color="blue", linewidth=line_width)
        plt.title("Tracking Error of $e_{\\beta}$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$e_{\\beta}$ / A',  fontdict=fontdict);


        # =================================
        #    FIG 3: Current (Alpha)
        # =================================
        plt.subplot(243)
        plt.plot(t, [x[0] for x in x_hist], label="$i_{\\alpha}$", color="blue", linewidth=line_width)
        plt.plot(t, [xd[0] for xd in r_hist], label="$i_{\\alpha}^*$", color="red", linewidth=line_width)
        plt.title("Tracking Result of $i_{\\alpha}$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$i_{\\alpha}$ / A',  fontdict=fontdict);

        # =================================
        #    FIG 4: Current (Beta)
        # =================================
        plt.subplot(244)
        plt.plot(t, [x[1] for x in x_hist], label="$i_{\\beta}$", color="blue", linewidth=line_width)
        plt.plot(t, [xd[1] for xd in r_hist], label="$i_{\\beta}^*$", color="red", linewidth=line_width)
        plt.title("Tracking Result of $i_{\\beta}$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$i_{\\beta}$ / A',  fontdict=fontdict);

        # =================================
        #    FIG 5: Voltage (Alpha)
        # =================================
        plt.subplot(245)
        plt.plot(t, [u[0] for u in u_hist], label="$v_\\alpha$", color="blue", linewidth=line_width)
        plt.title("Control Input", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$v_\\alpha$ / V',  fontdict=fontdict);
        
        # =================================
        #    FIG 5: Voltage (Alpha)
        # =================================
        plt.subplot(246)
        plt.plot(t, [u[1] for u in u_hist], label="$v_\\beta$", color="blue", linewidth=line_width)
        plt.title("Control Input", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$v_\\beta$ / V',  fontdict=fontdict);

        # =================================
        #    FIG 6: M Norm
        # =================================
        plt.subplot(247)
        plt.plot(t, M_hist, label="$||M||_2$", color="blue", linewidth=line_width)
        plt.title("Norm of $M$", fontdict=fontdict)
        plt.grid(True)
        plt.legend(fontsize=lgd_size)
        plt.xlabel('Time / s', fontdict=fontdict);
        plt.ylabel('$||M||_2$',  fontdict=fontdict);

        plt.show()

# ============================================
#         MAIN FUNCTION
# ============================================
if __name__ == "__main__":
    X_MAX = 10
    X_NUM = 1

    x_range = np.linspace(-X_MAX, X_MAX, X_NUM, dtype="f")

    _SAVE = {}

    for idx_x1 in range(len(x_range)):
        for idx_x2 in range(len(x_range)):
            x = np.array([
                x_range[idx_x1],
                x_range[idx_x2],
            ], dtype="f").reshape(2,1);
            xd0 = np.array([
                10.0, 
                -10.0, 
            ], dtype="f").reshape(2,1);

            print(f"Simulation {idx_x1}, {idx_x2}")
            print(f"x: {x.T}, xd0: {xd0.T}")
            x_hist, u_hist, r_hist, ud_hist, M_hist = simulate(x0=x, xd0=xd0)

            _SAVE[f'x{x}_xd{xd0}'] = {
                'x_hist': x_hist,
                'u_hist': u_hist,
                'r_hist': r_hist,
                'ud_hist': ud_hist,
                'M_hist': M_hist
            }

            if FIGURE_PLOT:
                result_plot()

    if RESULT_SAVE:
        with open('result.pickle', 'wb') as f:
            pickle.dump(_SAVE, f)
        print("Result saved to result.pickle")


        
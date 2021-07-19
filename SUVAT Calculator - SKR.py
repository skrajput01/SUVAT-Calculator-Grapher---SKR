# 05/07/21 - 19/07/21
from math import sqrt
import matplotlib.pyplot as plt
import numpy as np
""" PROGRAMMER'S NOTE: I am first and foremost a physicist - programming is a love and hobby of mine but the way my code is written reflects this. My calculations primarily use SUVAT rearrangements over iterative processes to obtain missing values. With SUVAT, there are some instances where two different outcomes are equally possible. The SUVAT calculator will plot both of these instances. The calculator also informs users of errors with input values such as negative time inputs, too few SUVAT inputs, etc. There were two goals I had when creating the calculator. The first was to be reliable 100% of the time, hence why there are so many if statements within if statements - which ensures correct calculations in all scenarios are performed, and when they are impossible to calculate, the user is informed of where the error stems from. The second goal was to create the most advanced calculator possible as I couldn't find anything quite as indepth as this one (which is able to calculate other important variables like max displacement and plot the potential trajectory of the projectile). 
"""
print("S U V A T Calculator - by Shyam Kumar Rajput", "\n")
print("PLEASE ONLY ENTER 3 KNOWN VARIABLES. LEAVE LINE BLANK FOR UNKNOWN VARIABLES", "\n")

# INITIAL INPUT VALUES [S, U, V] (Leave Line Blank if Unknown)
S = input("(Displacement)     S: ")
U = input("(Initial Velocity) U: ")
V = input("(Final Velocity)   V: ")
# INITIAL VALUES LETTER TEST [S, U, V] (Ensure Only Numerical Inputs)
test_S = S.isalpha()
test_U = U.isalpha()
test_V = V.isalpha()


# ERROR MESSAGES
def error(error_n):
    # ERROR 1: INVALID TIME CALCULATED
    if error_n == 1:
        print("ERROR[1]: Unable to Calculate a Valid Time. You May Have Entered an Incorrect Value")
    # ERROR 2: INVALID ACCELERATION CALCULATED
    elif error_n == 2:
        print("ERROR[2]: Unable to Calculate a Valid Acceleration. You May Have Entered an Incorrect Value")
    # ERROR 3: INVALID NON-NUMERICAL INPUT
    elif error_n == 3:
        print("ERROR[3]: You Have Entered a Letter for One or More Variables. Please Only Enter Numbers")
    # ERROR 4: INFINITE TIMES CALCULATED
    elif error_n == 4:
        print("ERROR[4]: An Infinite Number of Possible Times Were Calculated for These Variables")
    # ERROR 5: INVALID TIME INPUT
    elif error_n == 5:
        print("ERROR[5]: You Have Entered an Invalid Time")
    # ERROR 6: TOO FEW INPUTS
    elif error_n == 6:
        print("ERROR[6]: You Have Entered Too Few Variables. Please Enter 3 Known Variables")
    # ERROR 7: UNEXPECTED ERROR
    elif error_n == 7:
        print("ERROR[7]: An Unexpected Error Has Occurred")
    # ERROR 8: INVALID DISPLACEMENT CALCULATED
    elif error_n == 8:
        print("ERROR[8]: Unable to Calculate a Valid Displacement. You May Have Entered an Incorrect Value")
    # ERROR 9: INVALID FINAL VELOCITY CALCULATED
    elif error_n == 9:
        print("ERROR[9]: Unable to Calculate a Valid Final Velocity. You May Have Entered an Incorrect Value")
    # ERROR 10: INVALID Initial VELOCITY CALCULATED
    elif error_n == 10:
        print("ERROR[9]: Unable to Calculate a Valid Initial Velocity. You May Have Entered an Incorrect Value")
    elif error_n == 11:
        print("ERROR[11]: ADMIN ERROR")
    else:
        print("ADMIN ERROR")
    exit()


# KNOWN: S, U, V
if len(S) > 0 and len(U) > 0 and len(V) > 0:
    # TEST S, U, V: NUMERICAL INPUTS
    if test_S is False and test_U is False and test_V is False:
        print("\nFor Known: S, U, V")
        S = float(S)
        U = float(U)
        V = float(V)
        # CALCULATIONS
        # CALCULATE A
        try:
            A = round(((V**2) - (U**2))/(2*S), 2)
        except ZeroDivisionError:
            # ERROR: INVALID ACCELERATION CALCULATED
            error(2)
            A = None
        # CALCULATE T
        try:
            T = round((2*S)/(U + V), 2)
            # TEST T: VALID CALCULATIONS
            if T < 0:
                # ERROR: INVALID TIME CALCULATED
                error(1)
        except ZeroDivisionError:
            # ERROR: INVALID TIME CALCULATED
            error(1)
            T = None
    else:
        # ERROR: INVALID NON-NUMERICAL INPUTS
        error(3)
        T = None
        A = None
else:
    # SECONDARY INPUT VALUE [A]
    A = input("(Acceleration)     A: ")
    # SECONDARY VALUE LETTER TEST [A]
    test_A = A.isalpha()

    # KNOWN: S, U, A
    if len(S) > 0 and len(U) > 0 and len(A) > 0:
        # TEST S, U, A: NUMERICAL INPUTS
        if test_S is False and test_U is False and test_A is False:
            print("\nFor Known: S, U, A")
            S = float(S)
            U = float(U)
            A = float(A)
            # CALCULATIONS
            # TEST A: SUVAT SIMPLIFICATION (when A=0, S=U/T)
            if A == 0:
                try:
                    T = round(S/U, 2)
                    V = round(U, 2)
                    if T < 0:
                        # ERROR: INVALID TIME CALCULATED
                        error(1)
                except ZeroDivisionError:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                    T = None
            else:
                if (U**2 + 2*A*S) < 0:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                try:
                    T1 = (-U + sqrt((U**2) + (2*A*S)))/A
                    T2 = (-U - sqrt((U**2) + (2*A*S)))/A
                    # T1 and T2 are Unique and Positive
                    if (T1 >= 0 and T2 >= 0) and (T1 != T2):
                        print("In This Scenario, There are Two Possible Times & Therefore Two Possible Final Speeds")
                        V1 = round(U + A*T1, 2)
                        V2 = round(U + A*T2, 2)
                        T1 = round(T1, 2)
                        T2 = round(T2, 2)
                        print("S: ", S, "m")
                        print("U: ", U, "ms^-1")
                        print("V: ", V1, "ms^-1", "  (<- V1)")
                        print("   ", V2, "ms^-1", "  (<- V2)")
                        print("A: ", A, "ms^-2")
                        print("T: ", T1, "s", "      (<- T1)")
                        print("   ", T2, "s", "      (<- T2)\n")
                        # GRAPHS
                        time_1 = np.linspace(0, T1, 100)
                        time_2 = np.linspace(0, T2, 100)
                        # DISPLACEMENT-TIME GRAPH
                        # Set X-Axis & Y-Axis for Displacement
                        fig1, ax1 = plt.subplots()
                        ax1.spines["left"].set_position("zero")
                        ax1.spines["bottom"].set_position("zero")
                        # Displacement Plot
                        disp_1 = U*time_1 + 0.5*A*(time_1**2)
                        disp_2 = U*time_2 + 0.5*A*(time_2**2)
                        s_peak = round((-0.5*U**2)/A, 2)
                        t_peak = round(-U/A, 2)
                        t_comp = round((-2*U)/A, 2)
                        t_dashed = np.linspace(0, t_comp, 100)
                        s_dashed = U*t_dashed + 0.5*A*(t_dashed**2)
                        ax1.plot(time_1, disp_1, label=("S = {U}T + {hlf_A}T^2".format(U=U, hlf_A=0.5*A)), color="b")
                        ax1.plot(t_dashed, s_dashed, color="b", linestyle="dashed")
                        ax1.plot(time_2, disp_2, color="b")
                        ax1.scatter([0, T1, T2, t_peak, t_comp], [0, (U*T1 + 0.5*A*(T1**2)), (U*T2 + 0.5*A*(T2**2)), s_peak, 0], marker="x", color="k")
                        print("S at peak: ", s_peak, "m")
                        print("T at peak: ", t_peak, "s")
                        print("T to return:", t_comp, "s")
                        plt.title("Displacement-Time Graph")
                        plt.xlabel("Time [s]")
                        plt.ylabel("Displacement [m]")
                        plt.legend()
                        plt.grid()
                        plt.show()

                        # VELOCITY-TIME GRAPH
                        # Set X-Axis & Y-Axis for Velocity
                        fig2, ax2 = plt.subplots()
                        ax2.set_title("Velocity-Time Graph")
                        ax2.spines["top"].set_color("none")
                        ax2.spines["right"].set_color("none")
                        ax2.spines["left"].set_position("zero")
                        ax2.xaxis.set_ticks_position("bottom")
                        ax2.yaxis.set_ticks_position("left")
                        # Show Velocity Plot
                        vel_1 = np.linspace(U, V1, 100)
                        vel_2 = np.linspace(U, V2, 100)
                        ax2.plot(time_1, vel_1, label=("V = {A}T + {U}".format(A=A, U=U)), color="r")
                        ax2.plot(time_2, vel_2, color="r")
                        ax2.scatter([0, T1, T2], [U, U + A*T1, U + A*T2], marker="x", color="k")
                        if T1 > T2:
                            T = T1
                        else:
                            T = T2
                        plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
                        plt.xlabel("Time [s]")
                        plt.ylabel("Velocity [ms^-1]")
                        plt.legend()
                        ax2.grid()
                        plt.show()

                        # ACCELERATION-TIME GRAPH
                        # Set X-Axis & Y-Axis for Acceleration
                        fig3, ax3 = plt.subplots()
                        ax3.spines["top"].set_color("none")
                        ax3.spines["right"].set_color("none")
                        ax3.spines["left"].set_position("zero")
                        ax3.xaxis.set_ticks_position("bottom")
                        ax3.yaxis.set_ticks_position("left")
                        # Show Acceleration plots
                        if T1 > T2:
                            T = T1
                        else:
                            T = T2
                        time_acc = np.linspace(0, T, 100)
                        plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
                        acc = np.linspace(A, A, 100)
                        plt.plot(time_acc, acc, label=("A = {A}".format(A=A)), color="g")
                        ax3.scatter([0, T], [A, A], marker="x", color="k")
                        plt.title("Acceleration-Time Graph")
                        plt.xlabel("Time [s]")
                        plt.ylabel("Acceleration [ms^-2]")
                        plt.legend()
                        plt.grid()
                        plt.show()
                        exit()
                    # Only T1 is Positive
                    elif T1 >= 0:
                        V = round(U + A*T1, 2)
                        T = round(T1, 2)
                    # Only T2 is Positive
                    elif T2 >= 0:
                        V = round(U + A*T2, 2)
                        T = round(T2, 2)
                    # T1 and T2 are Negative
                    elif T1 < 0 and T2 < 0:
                        # ERROR: INVALID TIME CALCULATED
                        error(1)
                        T = None
                    else:
                        # ERROR: UNEXPECTED ERROR
                        error(7)
                        T = None
                except ZeroDivisionError:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                    T = None
        else:
            # ERROR: INVALID NON-NUMERICAL INPUT
            error(3)
            T = None

    # KNOWN: S, V, A
    elif len(S) > 0 and len(V) > 0 and len(A) > 0:
        # TEST S, V, A: NUMERICAL INPUTS
        if test_S is False and test_V is False and test_A is False:
            print("\nFor Known: S, V, A")
            S = float(S)
            V = float(V)
            A = float(A)
            # CALCULATIONS
            # TEST T, V: VALID T & V CALCULATIONS
            if (V**2 - 2*A*S) < 0:
                # ERROR: INVALID TIME CALCULATED
                error(1)
            # TEST A: SUVAT SIMPLIFICATION (when A=0, S=V/T)
            if A == 0:
                try:
                    T = round(S/V, 2)
                    U = round(V, 2)
                    if T < 0:
                        # ERROR: INVALID TIME CALCULATED
                        error(1)
                except ZeroDivisionError:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                    T = None
            else:
                if (V**2 - 2*A*S) < 0:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                try:
                    T1 = (V + sqrt((V**2) - (2*A*S)))/A
                    T2 = (V - sqrt((V**2) - (2*A*S)))/A
                    # T1 and T2 are Unique and Positive
                    if (T1 >= 0 and T2 >= 0) and (T1 != T2):
                        print("In This Scenario, There are Two Possible Times & Therefore Two Possible Final Speeds")
                        U1 = round(V - A*T1, 2)
                        U2 = round(V - A*T2, 2)
                        print("S: ", S, "m")
                        print("U: ", U1, "ms^-1", "  (<- U1)")
                        print("   ", U2, "ms^-1", "  (<- U2)")
                        print("V: ", V, "ms^-1")
                        print("A: ", A, "ms^-2")
                        print("T: ", round(T1, 2), "s", "      (<- T1)")
                        print("   ", round(T2, 2), "s", "      (<- T2)")
                        # GRAPHS
                        time_1 = np.linspace(0, T1, 100)
                        time_2 = np.linspace(0, T2, 100)
                        # DISPLACEMENT-TIME GRAPH
                        # Set X-Axis & Y-Axis for Displacement
                        fig4, ax4 = plt.subplots()
                        ax4.spines["left"].set_position("zero")
                        ax4.spines["bottom"].set_position("zero")
                        # Displacement Plot
                        disp_1 = V*time_1 - 0.5*A*(time_1**2)
                        disp_2 = V*time_2 - 0.5*A*(time_2**2)
                        s_peak = round((0.5*V**2)/A, 2)
                        t_peak = round(V/A, 2)
                        t_comp = round((2*V)/A, 2)
                        t_dashed = np.linspace(0, t_comp, 100)
                        s_dashed = V*t_dashed - 0.5*A*(t_dashed**2)
                        ax4.plot(time_1, disp_1, label=("S = {V}T + {mns_hlf_A}T^2".format(V=V, mns_hlf_A=-0.5*A)), color="b")
                        ax4.plot(time_2, disp_2, color="b")
                        ax4.plot(t_dashed, s_dashed, color="b", linestyle="dashed")
                        ax4.scatter([0, T1, T2, t_peak, t_comp], [0, (V*T1 - 0.5*A*(T1**2)), (V*T2 - 0.5*A*(T2**2)), s_peak, 0], marker="x", color="k")
                        print("S at peak: ", s_peak, "m")
                        print("T at peak: ", t_peak, "s")
                        print("T to return:", t_comp, "s")
                        plt.title("Displacement-Time Graph")
                        plt.xlabel("Time [s]")
                        plt.ylabel("Displacement [m]")
                        plt.legend()
                        plt.grid()
                        plt.show()

                        # VELOCITY-TIME GRAPH
                        # Set X-Axis & Y-Axis for Velocity
                        fig5, ax5 = plt.subplots()
                        ax5.set_title("Velocity-Time Graph")
                        ax5.spines["top"].set_color("none")
                        ax5.spines["right"].set_color("none")
                        ax5.spines["left"].set_position("zero")
                        ax5.xaxis.set_ticks_position("bottom")
                        ax5.yaxis.set_ticks_position("left")
                        # Show Velocity Plot
                        vel_1 = np.linspace(U1, V, 100)
                        vel_2 = np.linspace(U2, V, 100)
                        ax5.plot(time_1, vel_1, label=("U = {V} + {mns_A}T".format(mns_A=-A, V=V)), color="r")
                        ax5.plot(time_2, vel_2, color="r")
                        ax5.scatter([0, T1, T2], [U, V - A*T1, V - A*T2], marker="x", color="k")
                        if T1 > T2:
                            T = T1
                        else:
                            T = T2
                        plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
                        plt.xlabel("Time [s]")
                        plt.ylabel("Velocity [ms^-1]")
                        plt.legend()
                        ax5.grid()
                        plt.show()

                        # ACCELERATION-TIME GRAPH
                        # Set X-Axis & Y-Axis for Acceleration
                        fig6, ax6 = plt.subplots()
                        ax6.spines["top"].set_color("none")
                        ax6.spines["right"].set_color("none")
                        ax6.spines["left"].set_position("zero")
                        ax6.xaxis.set_ticks_position("bottom")
                        ax6.yaxis.set_ticks_position("left")
                        # Show Acceleration plots
                        if T1 > T2:
                            time_acc = np.linspace(0, T1, 100)
                            T = T1
                        else:
                            time_acc = np.linspace(0, T2, 100)
                            T = T2
                        plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
                        acc = np.linspace(A, A, 100)
                        plt.plot(time_acc, acc, label=("A = {A}".format(A=A)), color="g", marker="x")
                        ax6.scatter([0, T], [A, A], marker="x", color="k")
                        plt.title("Acceleration-Time Graph")
                        plt.xlabel("Time [s]")
                        plt.ylabel("Acceleration [ms^-2]")
                        plt.legend()
                        plt.grid()
                        plt.show()
                        exit()
                    # Only T1 is Positive
                    elif T1 >= 0:
                        U = round(V - A*T1, 2)
                        T = round(T1, 2)
                    # Only T2 is Positive
                    elif T2 >= 0:
                        U = round(V - A*T2, 2)
                        T = round(T2, 2)
                    # T1 and T2 are Negative
                    elif T1 < 0 and T2 < 0:
                        # ERROR: INVALID TIME CALCULATED
                        error(1)
                        T = None
                    else:
                        # ERROR: UNEXPECTED ERROR
                        error(7)
                        T = None
                except ZeroDivisionError:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
                    T = None
        else:
            # ERROR: INVALID NON-NUMERICAL INPUT
            error(3)
            T = None

    # KNOWN: U, V, A
    elif len(U) > 0 and len(V) > 0 and len(A) > 0:
        # TEST U, V, A: NUMERICAL INPUTS
        if test_U is False and test_V is False and test_A is False:
            print("\nFor Known: U, V, A")
            U = float(U)
            V = float(V)
            A = float(A)
            # CALCULATIONS
            # CALCULATE S
            try:
                S = round(((V**2) - (U**2))/(2*A), 2)
            except ZeroDivisionError:
                # ERROR: INVALID DISPLACEMENT CALCULATED
                error(8)
                T = None
            # CALCULATE T
            try:
                T = round((V - U)/A, 2)
                # TEST T: VALID CALCULATIONS
                if T < 0:
                    # ERROR: INVALID TIME CALCULATED
                    error(1)
            except ZeroDivisionError():
                # ERROR: INVALID TIME CALCULATED
                error(1)
                T = None
        else:
            # ERROR: INVALID NON-NUMERICAL INPUT
            error(3)
            T = None

    else:
        # THIRD INPUT VALUE [T]
        T = input("(Time)             T: ")
        # THIRD VALUE LETTER TEST [T]
        test_T = T.isalpha()

        # KNOWN: S, U, T
        if len(S) > 0 and len(U) > 0 and len(T) > 0:
            # TEST S, U, T: NUMERICAL INPUTS
            if test_S is False and test_U is False and test_T is False:
                print("\nFor Known: S, U, T")
                S = float(S)
                U = float(U)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATION
                # CALCULATE V
                try:
                    V = round(((2*S)/T) - U, 2)
                except ZeroDivisionError:
                    # ERROR: INVALID FINAL VELOCITY CALCULATED
                    error(9)
                    T = None
                # CALCULATE A
                try:
                    A = round(2*((S/(T**2)) - (U/T)), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID ACCELERATION CALCULATED
                    error(2)
                    T = None
            else:
                # ERROR: INVALID NON-NUMERICAL INPUT
                error(3)
                T = None

        # KNOWN: S, V, T
        elif len(S) > 0 and len(V) > 0 and len(T) > 0:
            # TEST S, V, T: NUMERICAL INPUTS
            if test_S is False and test_V is False and test_T is False:
                print("\nFor Known: S, V, T")
                S = float(S)
                V = float(V)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATION
                # CALCULATE U
                try:
                    U = round(((2*S)/T) - V, 2)
                except ZeroDivisionError:
                    # ERROR: INVALID INITIAL VELOCITY CALCULATED
                    error(10)
                    T = None
                # CALCULATE A
                try:
                    A = round(2 * ((V/T) - (S/(T**2))), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID ACCELERATION CALCULATED
                    error(2)
                    T = None
            else:
                # ERROR: INVALID NON-NUMERICAL INPUTS
                error(3)

        # KNOWN: S, A, T
        elif len(S) > 0 and len(A) > 0 and len(T) > 0:
            # TEST S, A, T: NUMERICAL INPUTS
            if test_S is False and test_A is False and test_T is False:
                print("\nFor Known: S, A, T")
                S = float(S)
                A = float(A)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATIONS
                # CALCULATE U
                try:
                    U = round(((S/T) - 0.5*A*T), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID INITIAL VELOCITY CALCULATED
                    error(10)
                    T = None
                # CALCULATE V
                try:
                    V = round(((S/T) + 0.5*A*T), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID FINAL VELOCITY CALCULATED
                    error(9)
                    T = None
            else:
                # ERROR: INVALID NON-NUMERICAL INPUTS
                error(3)

        # KNOWN: U, V, T
        elif len(U) > 0 and len(V) > 0 and len(T) > 0:
            # TEST U, V, T: NUMERICAL INPUTS
            if test_U is False and test_V is False and test_T is False:
                print("\nFor Known: U, V, T")
                U = float(U)
                V = float(V)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATION
                # CALCULATE S
                try:
                    S = round(0.5*(U + V)*T, 2)
                except ZeroDivisionError:
                    # ERROR: INVALID DISPLACEMENT CALCULATED
                    error(8)
                    T = None
                # CALCULATE A
                try:
                    A = round(((V - U)/T), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID ACCELERATION CALCULATED
                    error(2)
                    T = None
            else:
                # ERROR: INVALID NON-NUMERICAL INPUTS
                error(3)

        # KNOWN: U, A, T
        elif len(U) > 0 and len(A) > 0 and len(T) > 0:
            # TEST U, A, T: NUMERICAL INPUTS
            if test_U is False and test_A is False and test_T is False:
                print("\nFor Known: U, A, T")
                U = float(U)
                A = float(A)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATION
                # CALCULATE S
                try:
                    S = round((U*T + 0.5*A*(T**2)), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID DISPLACEMENT CALCULATED
                    error(8)
                    T = None
                # CALCULATE V
                try:
                    V = round((U + A*T), 2)
                except ZeroDivisionError:
                    # ERROR: INVALID FINAL VELOCITY CALCULATED
                    error(9)
                    T = None
            else:
                # ERROR: INVALID NON-NUMERICAL INPUTS
                error(3)

        # KNOWN: V, A, T
        elif len(V) > 0 and len(A) > 0 and len(T) > 0:
            # TEST V, A, T: NUMERICAL INPUTS
            if test_V is False and test_A is False and test_T is False:
                print("\nFor Known: V, A, T")
                V = float(V)
                A = float(A)
                T = float(T)
                # TEST T: VALID INPUT
                if T < 0:
                    # ERROR: INVALID TIME INPUT
                    error(5)
                # CALCULATION
                S = round((V*T - 0.5*A*(T**2)), 2)
                U = round((V - A*T), 2)
            else:
                # ERROR: INVALID NON-NUMERICAL INPUTS
                error(3)
        else:
            S = None
            U = None
            V = None
            A = None
            T = None
            # ERROR: TOO FEW INPUTS
            error(6)


# DISPLAY
print("S: ", S, "m")
print("U: ", U, "ms^-1")
print("V: ", V, "ms^-1")
print("A: ", A, "ms^-2")
print("T: ", T, "s\n")

# DISPLACEMENT-TIME GRAPH
time_disp = np.linspace(0, T, 100)
disp = U*time_disp + 0.5*A*(time_disp**2)
s_peak = round((-0.5*U**2)/A, 2)
t_peak = round(-U/A, 2)
t_comp = round(-2*U/A, 2)
t_dashed = np.linspace(0, t_comp, 100)
s_dashed = U*t_dashed + 0.5*A*(t_dashed**2)
fig7, ax7 = plt.subplots()
# Set X-Axis & Y-Axis for Displacement
ax7.spines["left"].set_position("zero")
ax7.spines["bottom"].set_position("zero")
# Show Displacement Plot
ax7.plot(time_disp, disp, label=("S = {U}T + {half_A}T^2".format(U=U, half_A=0.5*A)), color="b")
ax7.plot(t_dashed, s_dashed, color="b", linestyle="dashed")
ax7.scatter([0, t_peak, T, t_comp], [0, s_peak, S, 0], marker="x", color="k")
print("S at peak: ", s_peak, "m")
print("T at peak: ", t_peak, "s")
print("T to return:", t_comp, "s")
plt.title("Displacement-Time Graph")
plt.xlabel("Time [s]")
plt.ylabel("Displacement [m]")
plt.legend()
plt.grid()
plt.show()

# VELOCITY-TIME GRAPH
time_vel = np.linspace(0, T, 100)
vel = np.linspace(U, V, 100)
fig8, ax8 = plt.subplots()
# Set X-Axis & Y-Axis for Velocity
ax8.set_title("Velocity-Time Graph")
ax8.spines["top"].set_color("none")
ax8.spines["right"].set_color("none")
ax8.spines["left"].set_position("zero")
ax8.xaxis.set_ticks_position("bottom")
ax8.yaxis.set_ticks_position("left")
# Show Velocity Plot
ax8.plot(time_vel, vel, label=("V = {A}T + {U}".format(A=A, U=U)), color="r")
plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
ax8.scatter([0, T], [U, V], marker="x", color="k")
plt.legend()
plt.xlabel("Time [s]")
plt.ylabel("Velocity [ms^-1]")
ax8.grid()
plt.show()

# ACCELERATION-TIME GRAPH
time_acc = np.linspace(0, T, 100)
acc = np.linspace(A, A, 100)
fig9, ax9 = plt.subplots()
# Set X-Axis & Y-Axis for Acceleration
ax9.spines["top"].set_color("none")
ax9.spines["right"].set_color("none")
ax9.spines["left"].set_position("zero")
ax9.xaxis.set_ticks_position("bottom")
ax9.yaxis.set_ticks_position("left")
# Show Acceleration plots
plt.plot(time_acc, acc, label=("A = {A}".format(A=A)), color="g")
plt.hlines(y=0, xmin=0, xmax=T, color='k', linestyle='-')
ax9.scatter([0, T], [A, A], marker="x", color="k")
plt.title("Acceleration-Time Graph")
plt.xlabel("Time [s]")
plt.ylabel("Acceleration [ms^-2]")
plt.legend()
plt.grid()
plt.show()

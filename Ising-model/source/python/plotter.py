import matplotlib.pyplot as plt
import numpy as np

save_last_step = int(input("Do you want plot for each monte carlo cycles? yes = 1 and no = 0 " ))
crit = int(input("do you want to see critical point? yes = 1 and no = 0 ") )
#filename = str(input("Please enter a filename "))
file_array = np.genfromtxt("solution_lattice_n.txt")

faen = np.genfromtxt('solution_lattice_fuuck_dis_ordered.txt')
faen2 = np.genfromtxt('solution_lattice_fuuck2_dis_ordered.txt')

print(np.std(faen[:,0]))
print(np.var(faen2[:,0]))
#print(np.std(faen[:,0]))

#print(np.std(faen2[:,0]))

if save_last_step == True:
    file_array = []
    for i in range(4):
        file_array.append(np.genfromtxt("solution_lattices" + str(i) + ".txt"))

    def save_last_step_plotter(array_file):


        def Cv_analytical(T):
            k = 1
            J = 1
            beta = 1/(k*T)
            factor = 1/(k*T*T)

            a1 = 64*J*np.cosh(8*beta*J)/(np.cosh(8*J*beta) + 3)
            a2 = (8*J*np.sinh(8*beta*J)/(np.cosh(8*J*beta) + 3))**2

            return factor*(a1-a2)/4

        def X_analytical(T):
            k = 1
            J = 1
            beta = 1/(k*T)
            factor = 1/(k*T)
            a1 = (8*np.exp(8*beta*J) + 8)
            a2 = np.cosh(8*J*beta) + 3
            return factor*((a1/a2)/4)

        print(X_analytical(1))
        """
        temperature     = file_array[:,0]
        spin            = file_array[0,1]
        monte_carlo     = file_array[0,2]
        mean_E          = file_array[:,3]
        mean_E2         = file_array[:,4]
        Cv              = file_array[:,5]
        mean_M          = file_array[:,6]
        mean_M2         = file_array[:,7]
        X               = file_array[:,8]
        abs_M           = file_array[:,9]
        """
        #f = plt.figure(figsize=(10,3))

        plt.subplot(2,1,1)
        for i in range(0,4):
            plt.plot(file_array[i][:,0],file_array[i][:,5])
        plt.plot(file_array[0][:,0],Cv_analytical(file_array[0][:,0]))
        plt.title('Heat Capacity,with $L^{2} =4$')
        plt.ylabel('Heat capacity, $C_{v}$')
        plt.xlabel('temperature, $T$')
        plt.legend(['N = 5e4','N = 5e5', 'N = 5e6', 'analytical'], loc = 'best')
        plt.grid('on')
        plt.tight_layout()

        plt.subplot(2,1,2)
        for i in range(0,4):
            plt.plot(file_array[i][:,0],file_array[i][:,8])
        plt.plot(file_array[0][:,0],X_analytical(file_array[0][:,0]))
        plt.title('Magnetic Suscepbility, with $L^{2} =4$')
        plt.ylabel('Magnetic susceptbility $\\chi$')
        plt.xlabel('temperature, $T$')
        plt.legend(['N = 5e4','N = 5e5', 'N = 5e6', 'analytical'], loc = 'upper right')
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/analytical_vs_numerical_2x2_magsus_and_heat_cap.pdf')

        plt.show()

        plt.subplot(2,1,1)
        for i in range(0,4):
            plt.plot(file_array[i][:,0],file_array[i][:,3])
        plt.plot(file_array[0][:,0], 0.25*(-8*np.sinh(8/file_array[0][:,0]))/(np.cosh(8/file_array[0][:,0]) +3))
        plt.title("Expectation value ofenergy, with $L^{2} =4$")
        plt.ylabel('Mean energy $\\langle E \\rangle$')
        plt.xlabel('temperature $T$')
        plt.legend(['N = 5e4','N = 5e5', 'N = 5e6', 'analytical'], loc = 'best')
        plt.grid('on')
        plt.tight_layout()

        plt.subplot(2,1,2)
        for i in range(0,4):
            plt.plot(file_array[i][:,0],file_array[i][:,9])
        plt.plot(file_array[0][:,0], 0.25*(2*np.exp(8/file_array[0][:,0])+4)/(np.cosh(8/file_array[0][:,0]) +3))
        plt.title("Expectation value of magnetic moment, with $L^{2} =4$ ")
        plt.ylabel('$\\langle |M| \\rangle$ $[T]$')
        plt.xlabel('temperature $T$')
        plt.legend(['N = 5e4','N = 5e5', 'N = 5e6', 'analytical'], loc = 'best')
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/analytical_vs_numerical_2x2_mean_absmag_and_meanE.pdf')
        plt.show()



    save_last_step_plotter(file_array)
#else:

def save_each_plotter(ordered_state, disordered_state, arg1,arg2):
    plt.subplot(2,1,1)
    plt.semilogx(ordered_state[0][:,2],ordered_state[0][:,0], '-r')
    plt.semilogx(disordered_state[0][:,2],disordered_state[0][:,0], '--b')
    plt.title("Expectation value of energy, $L^{2} = $ " + arg1 + " lattice, for $T =$ " + arg2)
    plt.ylabel(' $\\langle E \\rangle$')
    plt.xlabel('Monte carlo cycles')
    plt.legend(['ordered state', 'disordered state'],loc="upper right")
    plt.grid('on')
    #plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/numerical_'+ arg1 +'_meanE' + arg2 + '.pdf')
    plt.tight_layout()
    plt.subplot(2,1,2)
    plt.semilogx(ordered_state[0][:,2],ordered_state[0][:,1], '-r')
    plt.semilogx(disordered_state[0][:,2],disordered_state[0][:,1], '--b')
    plt.title("Expectation value of Magnetic moment, $L^{2} = $ " + arg1 + " lattice, for $T =$ " + arg2)
    plt.ylabel(' $\\langle |M| \\rangle$')
    plt.xlabel('Monte carlo cycles')
    plt.legend(['ordered state', 'disordered state'], loc= "lower right")
    plt.grid('on')
    plt.tight_layout()
    #plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/numerical_'+ arg1 +'_mean_absM_and_meanE_' + arg2 + '_Temp.pdf')
    #plt.show()


    plt.subplot(2,1,1)
    plt.semilogx(ordered_state[0][:,2],ordered_state[0][:,3], '-r')
    plt.semilogx(disordered_state[0][:,2],disordered_state[0][:,3], '--b')
    plt.title("Accepted spin confiquration $\\xi$ , $L^{2}= $ " + arg1 + " lattice, for $T = 1.0$ ")
    plt.ylabel(' Accepted spin confiquration $\\xi$')
    plt.xlabel('Monte carlo cycles')
    plt.legend(['ordered state', 'disordered state'])
    plt.grid('on')
    plt.tight_layout()

    plt.subplot(2,1,2)
    plt.semilogx(ordered_state[1][:,2],ordered_state[1][:,3], '-r')
    plt.semilogx(disordered_state[1][:,2],disordered_state[1][:,3], '--b')
    plt.title("Accepted spin confiquration $\\xi$ , $L^{2}= $ " + arg1 + " lattice, for $T = 2.4$" )
    plt.ylabel(' Accepted spin confiquration $\\xi$')
    plt.xlabel('Monte carlo cycles')
    plt.legend(['ordered state', 'disordered state'])
    plt.grid('on')
    plt.tight_layout()
    plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/numerical_' + arg1 + '_flipps_E_and_M_subplots.pdf')
    plt.show()



def probability_distribution(ordered_state, disordered_state, arg1,arg2):
    energy_temp24, occourance24 = np.unique(ordered_state[1][:,-1], return_counts=True)
    energy_temp1, occourance1 = np.unique(ordered_state[0][:,-1], return_counts=True)

    plt.subplot(2,1,1)
    plt.title("Energy distribution, $ L^{2} = 400$ lattice, for T = 2.4")
    plt.xlabel("Energy")
    plt.ylabel("P(E)")
    prob24 = occourance24/np.sum(occourance24)
    plt.bar(energy_temp24,height = prob24, width = 0.01)
    plt.grid('on')
    plt.tight_layout()

    plt.subplot(2,1,2)

    plt.title("Energy distribution, $ L^{2} = 400$ lattice, for T = 1.0 ")
    plt.xlabel("Energy")
    plt.ylabel("P(E)")
    prob1 = occourance1/np.sum(occourance1)
    plt.bar(energy_temp1[0:4],height = prob1[0:4], width = 0.01)
    plt.grid('on')
    plt.tight_layout()
    plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/numerical_prob_dist_' + arg1 + '_flipps_subplot.pdf')
    plt.show()

if save_last_step == 2:
    disordered_state = []
    ordered_state = []
    lattice_size = str(input("lattice size LxL =  "))
    temperature = str (input("Choose temperature, T = " ))
    disordered_state.append(np.genfromtxt("solution_lattice1_dis_ordered.txt"))
    disordered_state.append(np.genfromtxt("solution_lattice_dis_ordered.txt"))
    ordered_state.append(np.genfromtxt("solution_lattice1_ordered.txt"))
    ordered_state.append(np.genfromtxt("solution_lattice_ordered.txt"))

    probability_distribution(ordered_state, disordered_state, lattice_size, temperature)
    save_each_plotter(ordered_state, disordered_state, lattice_size, temperature)

if crit == True:
    def crit_plott(list_file):
        tempCv = []
        maxvalCv = []

        tempX = []
        maxvalX = []
        for i in range(5):
            maxvalCv.append(np.argmax(list_file[i][:,5]))
            tempCv.append(list_file[i][maxvalCv[i],0])

            maxvalX.append(np.argmax(list_file[i][:,8]))
            tempX.append(list_file[i][maxvalX[i],0])
            print(tempX[i], list_file[i][maxvalX[i],5], list_file[i][maxvalX[i],8])


        plt.subplot(2,1,1)
        for i in range(5):
            plt.plot(list_file[i][:,0],list_file[i][:,5])

        plt.title('Heat Capacity for different lattice sizes with $T_{C}$')
        plt.ylabel('Heat capacity, $C_{v}$')
        plt.xlabel('temperature, $T$')
        plt.legend(['L = 40', 'L = 60', 'L = 80', 'L = 100', 'L = 120'])
        for i in range(5):
            plt.plot(tempCv[i],list_file[i][maxvalCv[i],5],'-o')
        #plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/analytical_vs_numerical_2x2_heatcap.pdf')
        plt.grid('on')
        plt.tight_layout()

        plt.subplot(2,1,2)
        for i in range(5):
            plt.plot(list_file[i][:,0],list_file[i][:,8])

        plt.title('Magnetic Suscepbility for different lattice sizes with $T_{C}$')
        plt.ylabel('Magnetic susceptbility $\\chi$')
        plt.xlabel('temperature, $T$')
        plt.legend(['L = 40', 'L = 60', 'L = 80', 'L = 100', 'L = 120'])
        for i in range(5):
            plt.plot(tempX[i],list_file[i][maxvalX[i],8], '-o')
        plt.grid('on')
        plt.tight_layout()
        plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/heatcap_magsus_subplot.pdf')
        plt.show()


        plt.subplot(2,1,1)
        for i in range(5):
            plt.plot(list_file[i][:,0],list_file[i][:,3])
        plt.title("Expectation value for energy for different lattice sizes with $T_{C}$")
        plt.ylabel('Mean energy $\\langle E \\rangle$')
        plt.xlabel('temperature $T$')
        plt.legend(['L = 40', 'L = 60', 'L = 80', 'L = 100', 'L = 120'])
        #plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/analytical_vs_numerical_2x2_meanE.pdf')
        plt.grid('on')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)

        plt.subplot(2,1,2)
        for i in range(5):
            plt.plot(list_file[i][:,0],list_file[i][:,9])
        plt.title("Expectation value for absolute magnetic moment\n for different lattice sizes with $T_{C}$")
        plt.ylabel('$\\langle |M| \\rangle$ $[T]$')
        plt.xlabel('temperature $T$')
        plt.legend(['L = 40', 'L = 60', 'L = 80', 'L = 100', 'L = 120'])
        plt.grid('on')
        plt.tight_layout(pad=0.6, w_pad=0.7, h_pad=1.5)
        plt.savefig('/Users/aramsalihi/Desktop/FYS4150/project4/figures/crit_subplot_meanE_meanM.pdf')
        plt.show()


    list_f = []
    list_f.append(np.loadtxt("solution_lattice_crit120.txt"))
    #list_f.append(np.loadtxt("solution_lattice_nnnn1.txt"))
    for i in range(40,130, 20):
        list_f.append(np.loadtxt("solution_lattice_crit" + str(i) + ".txt"))

    crit_plott(list_f)

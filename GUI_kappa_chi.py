#-*-coding:utf-8-*-
import numpy as np
# import numpy.linalg as LA
# import random
import tkinter
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg  import FigureCanvasTkAgg
from functools import partial
from tkinter import filedialog


def kappa1(chi, delta, p0, eta, expected_payoff):
    RE = expected_payoff[0]
    SE = expected_payoff[1]
    TE = expected_payoff[2]
    # PE = expected_payoff[3]
    mu = 1 - eta
    
    return RE - ((1 / delta - 1) * (1 - p0) + eta)/(mu - eta) * (TE - RE + (TE - SE) / (chi - 1))

def kappa2(chi, delta, p0, eta, expected_payoff):
    RE = expected_payoff[0]
    SE = expected_payoff[1]
    TE = expected_payoff[2]
    PE = expected_payoff[3]
    mu = 1 - eta
    
    return PE + eta / (mu - eta) * (PE - SE + (TE - SE) / (chi - 1)) \
            + (1 - delta) * p0 * (RE - PE + (TE - SE) / (chi - 1) + (mu * (TE - RE) - eta * (PE - SE)) / (mu - eta))

def kappa3(chi, delta, p0, eta, expected_payoff):
    RE = expected_payoff[0]
    SE = expected_payoff[1]
    TE = expected_payoff[2]
    PE = expected_payoff[3]
    mu = 1 - eta
    return RE - eta / (mu - eta) * (TE - RE + (TE - SE) / (chi - 1)) \
            - (1 - delta) * (1 - p0) * (RE - PE + (TE - SE) / (chi - 1) + (mu * (PE - SE) - eta * (TE - RE)) / (mu - eta))
                
def kappa4(chi, delta, p0, eta, expected_payoff):
    # RE = expected_payoff[0]
    SE = expected_payoff[1]
    TE = expected_payoff[2]
    PE = expected_payoff[3]
    mu = 1 - eta
    return PE + ((1 / delta - 1) * p0 +eta) / (mu - eta)* (PE - SE + (TE - SE) / (chi - 1))


def Quit():
    global root
    root.quit()
    root.destroy()

def change_5310(canvas, ax):
    global T,R,P,S
    if T == 5:
        T,R,P,S=1.5,1,0,-0.5
    else:
        T,R,P,S=5,3,1,0
    DrawCanvas(canvas, ax, colors = "gray")

def save_fig():
    filepath = filedialog.askdirectory(initialdir = dir)
    path = filepath + '\\fig.png'
    print('save_Image')
    print(path)
    plt.savefig(path)
    
def DrawCanvas(canvas, ax, colors = "gray"):
    ax.cla()
    
    delta = round(1-scale5.get()/100,2)
    p0 = round(scale_p0.get()/1000,3)
    eta = round(scale7.get()/1000,3)
    
    # plt.title(r"aaa", fontsize=15)
    RE = R*(1-eta)+S*eta
    SE = S*(1-eta)+R*eta
    TE = T*(1-eta)+P*eta
    PE = P*(1-eta)+T*eta
    mu = 1 -eta
    
    expected_payoff = [RE,SE,TE,PE]
    
    plt.ylabel(f"$\kappa$ ", fontsize=15)
    plt.xlabel(f"$\chi_c$ ", fontsize=15)
    
    plt.grid()
    
    plt.xlim([0, 20])
    plt.ylim([S, T])
        
    delta_c1 = (T - R) / (mu*(RE-PE)-eta*(TE-SE)+TE-RE)
    delta_c2 = (P - S) / (mu*(RE-PE)-eta*(TE-SE)+PE-SE)
    delta_c = max(delta_c1, delta_c2)
    
    #chi_c1 = 1 + (1 - delta) * (T - S) / (delta * (T - P) - (T - R))
    #chi_c2 = 1 + (1 - delta) * (T - S) / (delta * (R - S) - (P - S))
    chi_c1 = 1 + (1 - delta + 2 * delta * eta) * (TE - SE) / (delta * (mu*(RE-PE)-eta*(TE-SE)) - (1-delta)*(TE - RE))
    chi_c2 = 1 + (1 - delta + 2 * delta * eta) * (TE - SE) / (delta * (mu*(RE-PE)-eta*(TE-SE)) - (1-delta)*(PE - SE))
    chi_c = max(chi_c1, chi_c2)
    
    chi_list = [i for i in np.linspace(chi_c, 20, 1000) if i <= 20 and i > 1]
    kappa_list1 = []; kappa_list2 = []
    kappa_list3 = []; kappa_list4 = []
    kappa_p0_0_list = []; kappa_p0_1_list = []
    
    if delta >= delta_c:
        for chi in chi_list:
            kappa_list1.append(kappa1(chi, delta, p0, eta, expected_payoff))
            kappa_list2.append(kappa2(chi, delta, p0, eta, expected_payoff))
            
            kappa_list3.append(kappa3(chi, delta, p0, eta, expected_payoff))
            kappa_list4.append(kappa4(chi, delta, p0, eta, expected_payoff))
            
            
            
            kappa_p0_0_list.append(min(kappa1(chi, delta, 1, eta, expected_payoff),
                                       kappa3(chi, delta, 1, eta, expected_payoff)))
            kappa_p0_1_list.append(min(kappa2(chi, delta, 0, eta, expected_payoff),
                                       kappa4(chi, delta, 0, eta, expected_payoff)))
        
        plt.plot(chi_list, kappa_p0_0_list, 'k', markersize=3, alpha=0.4)
        plt.plot(chi_list, kappa_p0_1_list, 'k', markersize=3, alpha=0.4)
        
        plt.plot(chi_list, kappa_list1, 'g', markersize=3, label = 'kappa1')
        plt.plot(chi_list, kappa_list2, 'b', markersize=3, label = 'kappa2')
        plt.plot(chi_list, kappa_list3, 'r', markersize=3, label = 'kappa3')
        plt.plot(chi_list, kappa_list4, 'y', markersize=3, label = 'kappa4')
        
        
        plt.axvline(x=chi_c, markersize=3)
        plt.rcParams["font.size"] = 10
        plt.legend()
        
        txt.delete(0, tkinter.END)
        txt.insert(tkinter.END, r"exist ZD strategies")
                
        txt_Result.delete(0, tkinter.END)
        txt_Result.insert(tkinter.END, f'chi_c = {round(chi_c, 3)}')
        
        txt_Parameters.delete(0, tkinter.END)
        txt_Parameters.insert(tkinter.END, f'(T,R,P,S)=({T},{R},{P},{S}), p0={p0}, delta={delta}, eta={eta}')

    else:
        txt.delete(0, tkinter.END)
        txt.insert(tkinter.END, r"do not exist ZD strategies")
    
    canvas.draw()

if __name__ == "__main__":
    try:
        #generate GUI
        root = tkinter.Tk()
        root.geometry("700x550")
        root.title("GUI- kappa_chi_c")
        
        #generate graph
        fig,ax1 = plt.subplots(figsize = (8,6), dpi=70)
        # fig.gca().set_aspect('equal', adjustable='box')
        
        #generate Canvas
        Canvas = FigureCanvasTkAgg(fig, master=root)
        Canvas.get_tk_widget().grid(row=0, column=0, rowspan=1000)
        T,R,P,S = 1.5,1,0,-0.5
        # q_list=[[random.random(),random.random(),random.random(),random.random(),random.random()] for i in range(1000)]
        
        #ReDrawButton = tkinter.Button(text="Other Opponent", width=15, command=partial(change_q, Canvas, ax1))
        #ReDrawButton.grid(row=12, column=1, columnspan=1)
        TRPSButton = tkinter.Button(text="TRPS 5310", width=15, command=partial(change_5310, Canvas, ax1))
        TRPSButton.grid(row=15, column=1, columnspan=1)
        SaveButton = tkinter.Button(text="Save Fig", width=15, command=save_fig)
        SaveButton.grid(row=21, column=1, columnspan=1)
        # SaveButton = tkinter.Button(text="Method", width=15, command=partial(change_way_cal,Canvas, ax1))
        SaveButton.grid(row=10, column=1, columnspan=1)
        
        #scale0 = tkinter.Scale(root, label='p0', orient='h', from_=0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        #scale0.grid(row=2, column=3, columnspan=1)
        #scale1 = tkinter.Scale(root, label='p1', orient='h', from_=0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        #scale1.grid(row=3, column=3, columnspan=1)
        #scale2 = tkinter.Scale(root, label='p2',orient='h', from_=0.0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        #scale2.grid(row=4, column=3, columnspan=1)
        #scale3 = tkinter.Scale(root, label='p3',orient='h', from_=0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        #scale3.grid(row=5, column=3, columnspan=1)
        #scale4 = tkinter.Scale(root, label='p4',orient='h', from_=0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        #scale4.grid(row=6, column=3, columnspan=1)
        
        scale5 = tkinter.Scale(root, label='discout rate',orient='h', from_=0, to=100, command=partial(DrawCanvas, Canvas, ax1))
        scale5.grid(row=2, column=1, columnspan=1)
        scale_p0 = tkinter.Scale(root, label='p0',orient='h', from_=0, to=1000, command=partial(DrawCanvas, Canvas, ax1))
        scale_p0.grid(row=3, column=1, columnspan=1)
        scale7 = tkinter.Scale(root, label='eta', orient='h', from_=0, to=300, command=partial(DrawCanvas, Canvas, ax1))
        scale7.grid(row=4, column=1, columnspan=1)
        
        QuitButton = tkinter.Button(text="Quit", width=15, command=Quit)
        QuitButton.grid(row=30, column=1, columnspan=1)
        
        lbl_Parameters = tkinter.Label(text = 'Parameters:')
        lbl_Parameters.place(x=10, y=440)
        txt_Parameters = tkinter.Entry(width=50)
        txt_Parameters.place(x=10, y=460)
        # txt_Parameters.insert(tkinter.END, f'(T,R,P,S)=({T},{R},{P},{S})')
        
        lbl_Result = tkinter.Label(text = 'Result:')
        lbl_Result.place(x=350, y=440)
        txt_Result = tkinter.Entry(width=50)
        txt_Result.place(x=350, y=460)
        # txt_Result.insert(tkinter.END, f'(T,R,P,S)=({T},{R},{P},{S})')
        
        
        lbl = tkinter.Label(text = 'Message:')
        lbl.place(x=10, y=500)
        txt = tkinter.Entry(width=50)
        txt.place(x=10, y=520)
        txt.insert(tkinter.END, "Welcome!")
        
        DrawCanvas(Canvas,ax1)
        root.mainloop()
    except Exception as error:
        #import traceback
        #traceback.print_exc()
        txt.delete(0, tkinter.END)
        txt.insert(tkinter.END,error.__str__())
import matplotlib.pyplot as plt
import numpy as np
import os
import sys 
from EMS import *
from entangle import *

def picture(M, ir=[]): 
    if not M:
        print("The list M is empty. Nothing to plot.")
        return
        
    M = sort(M) 
    output_filename = "pic.jpg"

    A_coords = [item[0][0] for item in M]
    B_coords = [item[0][1] for item in M]

    if not A_coords or not B_coords:
        x_limit_max = 10.0
        y_limit_max = 10.0
        y_limit_min = 0.0
    else:
        max_A = max(A_coords)
        max_B = max(B_coords)
        min_B = min(B_coords) 
        padding = 1.5 
        x_limit_max = max(max_A + padding, 10.0)
        y_limit_max = max(max_B + padding, 5.0) 
        y_limit_min = min(0.0, min_B - padding) 

    fig, ax = plt.subplots(figsize=(10, 10))
    ax.set_title("The Graph associated to the Extended Multi-Segment", fontsize=14)
    ax.set_xlabel("A-coordinate (Positive Axis)", fontsize=12)
    ax.set_ylabel("B-coordinate", fontsize=12)

    ax.set_xlim(0, x_limit_max)
    ax.set_ylim(y_limit_min, y_limit_max)
    
    ax.set_xticks(np.arange(0, np.ceil(x_limit_max) + 1, 1))
    ax.set_yticks(np.arange(np.floor(y_limit_min), np.ceil(y_limit_max) + 1, 1))
    
    ax.grid(True, linestyle='--', alpha=0.6)
    
    # Draw Diagonals 
    x_diag = np.linspace(0, x_limit_max, 100) 
    ax.plot(x_diag, x_diag, 'k--', alpha=0.5, label='$y=x$')
    ax.plot(x_diag, -x_diag, 'k:', alpha=0.5, label='$y=-x$ (Extension)')
    # ----------------------------------------
    
    conv_M = conv2(M)
    A_halves = np.arange(0.5, x_limit_max + 1, 1.0)
    start_B = np.floor(y_limit_min) - 0.5
    B_halves = np.arange(start_B, y_limit_max + 1, 1.0)
    
    for A in A_halves:
        for B in B_halves:
            if B <= A and B >= -A:
                a = int(A + B + 1)
                b = int(A - B + 1)
                if Is_irred(0.5, (a, b), conv_M):
                    rect = plt.Rectangle((A - 0.5, B - 0.5), 1, 1, 
                                         facecolor='lightblue', alpha=0.5, 
                                         edgecolor='lightblue', linewidth=0.5, zorder=1)
                    ax.add_patch(rect)


    N = len(M)
    for i in range(N):
        for j in range(i + 1, N):

            (A_i, B_i), m_i = M[i]
            (A_j, B_j), m_j = M[j]
            if strict(M, i, j):
                diff_val = abs(A_i-A_j) + abs(B_i-B_j)- abs(diff(M, i, j))
                label_val = diff_val
                if label_val == 0:
                    rect_x = min(A_i, A_j)
                    rect_y = min(B_i, B_j)
                    rect_width = abs(A_i - A_j)
                    rect_height = abs(B_i - B_j)

                    rect = plt.Rectangle((rect_x, rect_y), rect_width, rect_height,
                                         facecolor='darkorchid', alpha=0.4, edgecolor='none', zorder=1)
                    ax.add_patch(rect)
                    print(f"Colored Purple Rectangle for points {i} and {j} (Diff={label_val})")
                ax.plot([A_i, A_j], [B_i, B_j], 'r-', alpha=0.6, zorder=2)
                mid_x = (A_i + A_j) / 2
                mid_y = (B_i + B_j) / 2
                ax.text(mid_x, mid_y, f'{label_val}', color='darkred', fontsize=10,
                        ha='center', va='center', bbox=dict(facecolor='white', alpha=0.7, edgecolor='none', boxstyle='round,pad=0.2'), zorder=3)
    for k, ((A, B), m) in enumerate(M):
        ax.plot(A, B, 'o', color='navy', markersize=8, zorder=4)
        ax.text(A + 0.1, B, f'{k+1}', fontsize=12, color='navy', ha='left', va='center', fontweight='bold')
    for [A,B] in ir:
        ax.plot(A, B, 'o', color='black', markersize=3, zorder=4)

    ax.set_aspect('equal', adjustable='box') 
    plt.savefig(output_filename)
    plt.close(fig) 
    print(f"\nPlot successfully saved to {output_filename}")

    try:
        if sys.platform.startswith('win'):
            os.startfile(output_filename)
        elif sys.platform.startswith('darwin'):
            os.system(f'open "{output_filename}"')
        elif sys.platform.startswith('linux'):
            os.system(f'xdg-open "{output_filename}"')
        else:
            print("Could not automatically open the file. Please open 'pic.jpg' manually.")
    except Exception as e:
        print(f"Error attempting to open file: {e}")
if __name__ == "__main__":
    picture(M)
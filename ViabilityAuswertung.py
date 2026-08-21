import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import filedialog
import os

def parse_lod(val):
    """Parses strings like '10^4' into numeric floats (10000.0)."""
    if isinstance(val, str) and '^' in val:
        base, exp = val.split('^')
        return float(base) ** float(exp)
    return float(val)

def main():
    # 1. Dateipfad des Skripts ermitteln und Subfolder für den Output erstellen
    try:
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        script_dir = os.getcwd()
        
    output_dir = os.path.join(script_dir, 'Plots_Output')
    os.makedirs(output_dir, exist_ok=True)

    # 2. Explorer-Fenster zur Dateiauswahl öffnen
    root = tk.Tk()
    root.withdraw()
    file_path = filedialog.askopenfilename(
        title="Select your Excel file",
        filetypes=[("Excel files", "*.xlsx *.xls")]
    )
    
    if not file_path:
        print("No file selected. Exiting.")
        return

    # 3. Daten anhand der Spaltenindizes einlesen (0-basiert)
    # B(1): Time, C(2): CFU_M1, D(3): LOD_M1, E(4): CFU_FeCN, F(5): LOD_FeCN, G(6): CFU_ABTS, H(7): LOD_ABTS
    # I(8): FeCN_OD600, L(11): OD420_corr, M(12): M1_OD600, R(17): ABTS_OD600
    df = pd.read_excel(file_path, usecols=[1, 2, 3, 4, 5, 6, 7, 8, 11, 12, 17])
    
    df.columns = [
        'Time', 'CFU_M1', 'LOD_M1', 'CFU_FeCN', 'LOD_FeCN', 
        'CFU_ABTS', 'LOD_ABTS', 'OD600_FeCN', 'OD420_corr', 'OD600_M1', 'OD600_ABTS'
    ]
    df['Time'] = df['Time'] / 24.0

    # 4. LOD-Werte formatieren
    for col in ['LOD_M1', 'LOD_FeCN', 'LOD_ABTS']:
        df[col] = df[col].apply(parse_lod)

    # 5. 0-Werte auf die jeweilige LOD setzen
    df['CFU_M1_plot'] = np.where(df['CFU_M1'] == 0, df['LOD_M1'], df['CFU_M1'])
    df['CFU_FeCN_plot'] = np.where(df['CFU_FeCN'] == 0, df['LOD_FeCN'], df['CFU_FeCN'])
    df['CFU_ABTS_plot'] = np.where(df['CFU_ABTS'] == 0, df['LOD_ABTS'], df['CFU_ABTS'])

    # 6. Mittelwerte und Standardabweichungen berechnen
    df_mean = df.groupby('Time').mean().reset_index()
    df_std = df.groupby('Time').std().reset_index()

    # --- Plot-Einstellungen ---
    colors = {'M1': '#0072B2', 'FeCN': '#E69F00', 'ABTS': '#CC79A7'}
    labels = {'M1': 'M1', 'FeCN': '5 mM FeCN', 'ABTS': '100 \u03BCM ABTS'}
    
    # Hilfsfunktion für CFU-Daten
    def plot_cfu_data(ax, key):
        y = df_mean[f'CFU_{key}_plot']
        yerr_upper = df_std[f'CFU_{key}_plot']
        
        # DER WISSENSCHAFTLICHE FIX:
        max_lower_err = y - df_mean[f'LOD_{key}']
        yerr_lower = np.clip(yerr_upper, 0, max_lower_err)
        
        ax.errorbar(
            df_mean['Time'], y, yerr=[yerr_lower, yerr_upper],
            fmt='-o', color=colors[key], capsize=4, label=labels[key], zorder=3
        )
        
        # LOD-Indikator (Dreiecke nach unten)
        lod_mask = df_mean[f'CFU_{key}'] == 0
        if lod_mask.any():
            ax.plot(
                df_mean.loc[lod_mask, 'Time'], df_mean.loc[lod_mask, f'CFU_{key}_plot'],
                linestyle='None', marker='v', color=colors[key], markersize=9, 
                markeredgecolor='black', markerfacecolor='white', zorder=5
            )

    # Hilfsfunktion für OD600-Daten 
    def plot_od600_data(ax, key):
        ax.errorbar(
            df_mean['Time'], df_mean[f'OD600_{key}'], yerr=df_std[f'OD600_{key}'],
            fmt='-x', color=colors[key], capsize=4, label=labels[key] 
        )

    # Dummy-Eintrag für die Legende (LOD Indikator)
    lod_dummy = plt.Line2D([0], [0], marker='v', color='w', markeredgecolor='black', markerfacecolor='white', markersize=9)

    # ==========================================
    # Bild 1: CFU (log) & OD420 korrigiert
    # ==========================================
    fig1, ax1 = plt.subplots(figsize=(10, 6))
    ax2 = ax1.twinx()

    plot_cfu_data(ax1, 'M1')
    plot_cfu_data(ax1, 'FeCN')
    plot_cfu_data(ax1, 'ABTS')

    ax2.errorbar(
        df_mean['Time'], df_mean['OD420_corr'], yerr=df_std['OD420_corr'],
        fmt='--D', color='black', capsize=4, label='FeCN (OD$_{420}$ corrected)'
    )

    ax1.set_yscale('log')
    ax1.set_ylim(bottom=1) 
    ax1.set_xlabel('Time (days)', fontsize=12)
    ax1.set_ylabel('CFU/mL', fontsize=12, color='black')
    
    # Y-Achsen-Label für OD420 angepasst
    ax2.set_ylabel('FeCN (OD$_{420}$ corrected)', fontsize=12, color='black')
    
    lines1, labels_leg1 = ax1.get_legend_handles_labels()
    lines2, labels_leg2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2 + [lod_dummy], labels_leg1 + labels_leg2 + ['≤ Limit of Detection (LOD)'], loc='upper right')
    fig1.tight_layout()

    # ==========================================
    # Bild 2: Nur CFU (log)
    # ==========================================
    fig2, ax3 = plt.subplots(figsize=(10, 6))

    plot_cfu_data(ax3, 'M1')
    plot_cfu_data(ax3, 'FeCN')
    plot_cfu_data(ax3, 'ABTS')

    ax3.set_yscale('log')
    ax3.set_ylim(bottom=1)
    ax3.set_xlabel('Time (days)', fontsize=12)
    ax3.set_ylabel('CFU/mL', fontsize=12)
    ax3.legend(lines1 + [lod_dummy], labels_leg1 + ['≤ Limit of Detection (LOD)'], loc='upper right')
    fig2.tight_layout()

    # ==========================================
    # Bild 3: Nur OD420 korrigiert
    # ==========================================
    fig3, ax4 = plt.subplots(figsize=(10, 6))
    
    ax4.errorbar(
        df_mean['Time'], df_mean['OD420_corr'], yerr=df_std['OD420_corr'],
        fmt='-D', color='black', capsize=4, label='FeCN (OD$_{420}$ corrected)'
    )

    ax4.set_xlabel('Time (days)', fontsize=12)
    
    # Y-Achsen-Label für OD420 angepasst
    ax4.set_ylabel('FeCN (OD$_{420}$ corrected)', fontsize=12)
    ax4.legend(loc='upper right')
    fig3.tight_layout()

    # ==========================================
    # Bild 4: Nur OD600 (Zellwachstum / Trübung)
    # ==========================================
    fig4, ax5 = plt.subplots(figsize=(10, 6))
    
    plot_od600_data(ax5, 'M1')
    plot_od600_data(ax5, 'FeCN')
    plot_od600_data(ax5, 'ABTS')

    ax5.set_xlabel('Time (days)', fontsize=12)
    ax5.set_ylabel('Optical Density (OD$_{600}$)', fontsize=12)
    ax5.legend(loc='upper right')
    fig4.tight_layout()

    # --- Bilder speichern ---
    path_fig1 = os.path.join(output_dir, 'Fig1_CFU_and_OD420.png')
    path_fig2 = os.path.join(output_dir, 'Fig2_CFU_only.png')
    path_fig3 = os.path.join(output_dir, 'Fig3_OD420_only.png')
    path_fig4 = os.path.join(output_dir, 'Fig4_OD600_only.png')
    
    fig1.savefig(path_fig1, dpi=600, bbox_inches='tight')
    fig2.savefig(path_fig2, dpi=600, bbox_inches='tight')
    fig3.savefig(path_fig3, dpi=600, bbox_inches='tight')
    fig4.savefig(path_fig4, dpi=600, bbox_inches='tight')
    
    print(f"Abbildungen wurden erfolgreich gespeichert in:\n{output_dir}")

    # --- Zeige Plots auf dem Bildschirm ---
    plt.show()

if __name__ == "__main__":
    main()
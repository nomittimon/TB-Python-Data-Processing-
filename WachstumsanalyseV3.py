import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
from scipy import stats
import os
import csv
import itertools
import datetime
import sys
import tkinter as tk
from tkinter import filedialog, ttk, messagebox

# ==========================================
# 1. INTERAKTIVE DATEIAUSWAHL & SETTINGS
# ==========================================
root = tk.Tk()
root.withdraw()

print("Bitte wähle die Optical Density Datei aus...")
input_datei = filedialog.askopenfilename(
    title="Wähle die Optical Density Datei aus", 
    filetypes=[("CSV Dateien", "*.csv"), ("Alle Dateien", "*.*")]
)

if not input_datei:
    print("Keine Datei ausgewählt. Abbruch.")
    sys.exit()

output_ordner = "Ergebnisse_Diagramme"
os.makedirs(output_ordner, exist_ok=True)

try:
    zeitstempel = os.path.getmtime(input_datei)
    datei_datum = datetime.datetime.fromtimestamp(zeitstempel).strftime('%y%m%d') + "_"
except Exception as e:
    print(f"Hinweis: Metadaten konnten nicht gelesen werden ({e}). Verwende leeren Präfix.")
    datei_datum = ""

# ==========================================
# 2. READ PLATE MAP & IDENTIFY GROUPS
# ==========================================
gruppen_zuordnung = {}

with open(input_datei, newline='', encoding='utf-8') as fh:
    rows = list(csv.reader(fh))

start_idx = None
for i, r in enumerate(rows):
    if r and any(cell.strip() == 'Group Label' for cell in r):
        start_idx = i
        break

if start_idx is None:
    raise ValueError("'Group Label' wurde in der Eingabedatei nicht gefunden.")

group_labels = rows[start_idx]
well_names   = rows[start_idx + 1]

mediatoren = ['ABTS', 'PMS', 'PES', 'FeCN']

for label, well in zip(group_labels, well_names):
    label, well = label.strip(), well.strip()

    if not label or label in ('Group Label', 'Group Averages', 'Group Standard Deviations'):
        continue
    if len(well) < 2 or well[0] not in "ABCDEFGH" or not well[1:].isdigit():
        continue 

    if label == 'Vnat_WT':
        label = 'Vnat_dns'
        
    if label == 'Vnat_MtrB':
        label = 'Vnat_MtrAB'

    for med in mediatoren:
        if label == f'Vnat_{med}_MtrB':
            label = f'Vnat_{med}_MtrAB'

    if label not in gruppen_zuordnung:
        gruppen_zuordnung[label] = []
    gruppen_zuordnung[label].append(well)

blank_gruppen = (
    'M1_Blank', 'ABTS_M1_Blank', 'PMS_M1_Blank', 'PES_M1_Blank', 'FeCN_M1_Blank',
    'Vnat_Blank', 'Vnat_ABTS_Blank', 'Vnat_PMS_Blank', 'Vnat_PES_Blank', 'Vnat_FeCN_Blank'
)
experiment_gruppen = [g for g in gruppen_zuordnung.keys() if g not in blank_gruppen]

root_select = tk.Toplevel(root)
root_select.title("Gruppen aus der Analyse entfernen")
root_select.geometry("450x520")

lbl = tk.Label(root_select, text="Das Programm hat folgende Gruppen identifiziert.\nWähle alle Gruppen aus, die GESTRICHEN werden sollen:", justify="left", font=("Arial", 10, "bold"))
lbl.pack(pady=15, padx=15)

frame = tk.Frame(root_select)
frame.pack(pady=5, fill="both", expand=True, padx=20)

scrollbar = tk.Scrollbar(frame, orient="vertical")
listbox = tk.Listbox(frame, selectmode=tk.MULTIPLE, yscrollcommand=scrollbar.set, font=("Arial", 10))
scrollbar.config(command=listbox.yview)

scrollbar.pack(side="right", fill="y")
listbox.pack(side="left", fill="both", expand=True)

for grp in sorted(experiment_gruppen):
    listbox.insert(tk.END, grp)

gestrichene_gruppen = []

def bestaetigen():
    global gestrichene_gruppen
    gewaehlte_indizes = listbox.curselection()
    gestrichene_gruppen = [listbox.get(i) for i in gewaehlte_indizes]
    root_select.destroy()

btn = tk.Button(root_select, text="Auswahl bestätigen", command=bestaetigen, font=("Arial", 10, "bold"), bg="#d3d3d3", fg="black")
btn.pack(pady=20)

root.wait_window(root_select)

if gestrichene_gruppen:
    for grp in gestrichene_gruppen:
        if grp in gruppen_zuordnung:
            del gruppen_zuordnung[grp]

display_names = {}
verbleibende_gruppen = [g for g in gruppen_zuordnung.keys() if g not in blank_gruppen]

if verbleibende_gruppen:
    root_rename = tk.Toplevel(root)
    root_rename.title("Gruppen umbenennen")
    root_rename.geometry("750x500")

    lbl_ren = tk.Label(root_rename, text="Benenne die Gruppen für die Diagramme um:", font=("Arial", 10, "bold"))
    lbl_ren.pack(pady=15, padx=10)

    frame_ren = tk.Frame(root_rename)
    frame_ren.pack(fill="both", expand=True, padx=15)

    canvas = tk.Canvas(frame_ren)
    scrollbar_ren = tk.Scrollbar(frame_ren, orient="vertical", command=canvas.yview)
    scrollable_frame = tk.Frame(canvas)

    scrollable_frame.bind("<Configure>", lambda e: canvas.configure(scrollregion=canvas.bbox("all")))
    canvas.create_window((0, 0), window=scrollable_frame, anchor="nw")
    canvas.configure(yscrollcommand=scrollbar_ren.set)

    canvas.pack(side="left", fill="both", expand=True)
    scrollbar_ren.pack(side="right", fill="y")

    entries = {}
    for idx, grp in enumerate(sorted(verbleibende_gruppen)):
        tk.Label(scrollable_frame, text=grp, font=("Arial", 10)).grid(row=idx, column=0, padx=10, pady=5, sticky="w")
        ent = tk.Entry(scrollable_frame, width=55, font=("Arial", 10))
        
        default_text = grp
        if grp.startswith("Vnat_"):
            parts = grp.split("_")
            if "WT" in grp: 
                rest = "_".join(parts[1:])
                default_text = rf"$\mathit{{V. natriegens}}$ {rest}"
            else:
                if len(parts) >= 3:
                    mediator = parts[1]
                    gen = "_".join(parts[2:])
                    gen = gen[0].lower() + gen[1:]
                    default_text = rf"$\mathit{{V. natriegens}}$ $\mathit{{\Delta {gen}}}$ {mediator}"
                elif len(parts) == 2:
                    gen = parts[1]
                    default_text = rf"$\mathit{{V. natriegens}}$ $\mathit{{\Delta {gen}}}$"
                
        ent.insert(0, default_text)
        ent.grid(row=idx, column=1, padx=10, pady=5)
        entries[grp] = ent

    def speichern_umbenennung():
        for grp, ent in entries.items():
            display_names[grp] = ent.get()
        root_rename.destroy()

    btn_ren = tk.Button(root_rename, text="Umbenennung bestätigen & Skript ausführen", command=speichern_umbenennung, font=("Arial", 10, "bold"), bg="#d3d3d3")
    btn_ren.pack(side="bottom", pady=20)

    root.wait_window(root_rename)

# ==========================================
# 3. READ DATA & BLANK CORRECTION
# ==========================================
df = pd.read_csv(input_datei, skiprows=start_idx + 1)
zeit = pd.to_numeric(df['Duration (Hours)'], errors='coerce')

wells_m1_blank   = gruppen_zuordnung.get('Vnat_Blank', gruppen_zuordnung.get('M1_Blank', []))
wells_abts_blank = gruppen_zuordnung.get('Vnat_ABTS_Blank', gruppen_zuordnung.get('ABTS_M1_Blank', []))
wells_pms_blank  = gruppen_zuordnung.get('Vnat_PMS_Blank', gruppen_zuordnung.get('PMS_M1_Blank', []))
wells_pes_blank  = gruppen_zuordnung.get('Vnat_PES_Blank', gruppen_zuordnung.get('PES_M1_Blank', []))
wells_fecn_blank = gruppen_zuordnung.get('Vnat_FeCN_Blank', gruppen_zuordnung.get('FeCN_M1_Blank', []))

mean_m1_blank   = df[wells_m1_blank].mean(axis=1) if wells_m1_blank else 0
mean_abts_blank = df[wells_abts_blank].mean(axis=1) if wells_abts_blank else 0
mean_pms_blank  = df[wells_pms_blank].mean(axis=1) if wells_pms_blank else 0
mean_pes_blank  = df[wells_pes_blank].mean(axis=1) if wells_pes_blank else 0
mean_fecn_blank = df[wells_fecn_blank].mean(axis=1) if wells_fecn_blank else 0

for grp, wells in gruppen_zuordnung.items():
    if grp in blank_gruppen:
        continue
    if 'ABTS' in grp: baseline = mean_abts_blank
    elif 'PMS' in grp: baseline = mean_pms_blank
    elif 'PES' in grp: baseline = mean_pes_blank
    elif 'FeCN' in grp: baseline = mean_fecn_blank
    else: baseline = mean_m1_blank
    df[wells] = df[wells].sub(baseline, axis=0)

# ==========================================
# 4. KINETISCHE WACHSTUMSPARAMETER (EMPFEHLUNG)
# ==========================================
def berechne_kinetik_empfehlung(zeit_series, od_series, window_size=7): 
    od_smooth = od_series.rolling(window=3, center=True, min_periods=1).mean()
    max_od = od_smooth.max()
    
    basis_od_linear = od_smooth.iloc[:5].mean()
    basis_od_linear = max(0.005, basis_od_linear) 
    y_0_log = np.log(basis_od_linear)
    
    untere_grenze = max(0.05, basis_od_linear + 0.04) 
    obere_grenze = max_od * 0.75 
    
    mask = (od_smooth >= untere_grenze) & (od_smooth <= obere_grenze)
    
    if mask.sum() < window_size:
        mask = od_smooth > (basis_od_linear + 0.02)
        if mask.sum() < window_size:
            return None, None, None, None, None, None
            
    t = zeit_series[mask].values
    y = od_smooth[mask].values
    ln_y = np.log(y)
    
    mu_max = 0
    best_intercept = 0
    best_t_window = None
    
    for i in range(len(t) - window_size + 1):
        t_window = t[i:i+window_size]
        ln_y_window = ln_y[i:i+window_size]
        
        slope, intercept, r_value, _, _ = stats.linregress(t_window, ln_y_window)
        
        if slope > mu_max and r_value**2 > 0.95:
            if y[i+window_size-1] > y[i] * 1.20:
                mu_max = slope
                best_intercept = intercept
                best_t_window = t_window
            
    if mu_max == 0:
        for i in range(len(t) - window_size + 1):
            t_window = t[i:i+window_size]
            ln_y_window = ln_y[i:i+window_size]
            slope, intercept, r_value, _, _ = stats.linregress(t_window, ln_y_window)
            if slope > mu_max and r_value**2 > 0.85:
                if y[i+window_size-1] > y[i] * 1.15:
                    mu_max = slope
                    best_intercept = intercept
                    best_t_window = t_window

    if mu_max <= 0:
        return None, None, None, None, None, None
        
    t_d = np.log(2) / mu_max
    lag_time = max(0, (y_0_log - best_intercept) / mu_max)
    
    return lag_time, t_d, mu_max, best_intercept, best_t_window, y_0_log

# ==========================================
# 4.5. INTERAKTIVE KINETIK KLASSE (UI)
# ==========================================
class InteraktiverFitter:
    def __init__(self, t, y, name, rec_lag, rec_td, rec_mu, rec_b, rec_t_window, rec_y0):
        self.t = t
        self.y = y
        self.name = name
        
        self.y_smooth = pd.Series(self.y).rolling(window=3, center=True, min_periods=1).mean().values
        self.ln_y_smooth = np.log(np.maximum(self.y_smooth, 1e-6))
        
        self.final_lag = rec_lag
        self.final_td = rec_td
        self.current_mu = rec_mu
        self.current_b = rec_b
        
        if rec_y0 is not None:
            self.baseline_ln = rec_y0
        else:
            self.baseline_ln = np.log(max(0.005, np.mean(self.y_smooth[:5]))) if len(self.y_smooth)>=5 else -5.3

        self.fig, self.ax = plt.subplots(figsize=(10, 7))
        plt.subplots_adjust(bottom=0.35)
        self.fig.canvas.manager.set_window_title(f"Interaktive Kinetik: {self.name}")
        
        self.ax.plot(self.t, np.log(np.maximum(self.y, 1e-6)), 'ko', markersize=3, alpha=0.3, label="Rohdaten ln(OD)")
        self.ax.plot(self.t, self.ln_y_smooth, 'k-', linewidth=1.5, alpha=0.8, label="Geglättet ln(OD)")
        
        self.exp_patch = None
        self.base_patch = None
        
        if rec_t_window is not None:
            self.exp_patch = self.ax.axvspan(rec_t_window[0], rec_t_window[-1], color='green', alpha=0.2, label='Empf. Exp-Phase')
        
        self.ax.set_title(f"{self.name}\n1. Bereich markieren  2. Zuweisen ('Exp-Phase' oder 'Baseline')")
        self.ax.set_xlabel("Zeit [h]")
        self.ax.set_ylabel("ln(OD)")
        
        self.selected_t = None
        self.selected_y = None
        
        self.line_exp, = self.ax.plot([], [], 'r-', lw=2, label="Aktueller Exp-Fit")
        self.line_base = self.ax.axhline(self.baseline_ln, color='blue', linestyle='--', alpha=0.5, label="Aktuelle Baseline")
        self.ax.legend(loc="upper left")
        
        self.span = widgets.SpanSelector(
            self.ax, self.on_select, 'horizontal', useblit=True, 
            interactive=True, props=dict(alpha=0.3, facecolor='gray')
        )
        
        ax_btn_exp = plt.axes([0.1, 0.15, 0.25, 0.075])
        self.btn_exp = widgets.Button(ax_btn_exp, 'Markierung als Exp-Phase')
        self.btn_exp.on_clicked(self.set_exp)
        
        ax_btn_base = plt.axes([0.4, 0.15, 0.25, 0.075])
        self.btn_base = widgets.Button(ax_btn_base, 'Markierung als Baseline (Lag)')
        self.btn_base.on_clicked(self.set_base)
        
        ax_btn_done = plt.axes([0.7, 0.15, 0.2, 0.075])
        self.btn_done = widgets.Button(ax_btn_done, 'Bestätigen & Weiter', color='#d3d3d3')
        self.btn_done.on_clicked(self.done)
        
        self.info_text = self.fig.text(0.1, 0.05, "", fontsize=12, fontweight='bold')
        self.update_plot_with_fit()
        
    def on_select(self, xmin, xmax):
        mask = (self.t >= xmin) & (self.t <= xmax)
        self.selected_t = self.t[mask]
        self.selected_y = self.ln_y_smooth[mask] 
        
    def set_exp(self, event):
        if self.selected_t is not None and len(self.selected_t) > 1:
            slope, intercept, _, _, _ = stats.linregress(self.selected_t, self.selected_y)
            if slope > 0:
                self.current_mu = slope
                self.current_b = intercept
                
                if self.exp_patch is not None:
                    self.exp_patch.remove()
                self.exp_patch = self.ax.axvspan(self.selected_t[0], self.selected_t[-1], color='green', alpha=0.2)
                
                self.update_plot_with_fit()
                
    def set_base(self, event):
        if self.selected_t is not None and len(self.selected_t) > 0:
            self.baseline_ln = np.mean(self.selected_y)
            
            if self.base_patch is not None:
                self.base_patch.remove()
            self.base_patch = self.ax.axvspan(self.selected_t[0], self.selected_t[-1], color='blue', alpha=0.15)
            
            self.line_base.set_ydata([self.baseline_ln, self.baseline_ln])
            self.update_plot_with_fit()
            
    def update_plot_with_fit(self):
        if self.current_mu is not None and self.current_mu > 0:
            y_fit = self.current_mu * self.t + self.current_b
            self.line_exp.set_data(self.t, y_fit)
            self.final_td = np.log(2) / self.current_mu
            self.final_lag = max(0, (self.baseline_ln - self.current_b) / self.current_mu)
        else:
            self.line_exp.set_data([], [])
            self.final_td = None
            self.final_lag = None
            
        td_str = f"{self.final_td:.2f} h" if self.final_td else "N/A"
        lag_str = f"{self.final_lag:.2f} h" if self.final_lag else "N/A"
        self.info_text.set_text(f"Aktuelle Werte -> Verdopplungszeit: {td_str} | Lag-Phase: {lag_str}")
        self.fig.canvas.draw_idle()
        
    def done(self, event):
        plt.close(self.fig)
        
    def run(self):
        plt.show(block=True)
        return self.final_lag, self.final_td

# ==========================================
# 5. PLOTTING FUNKTION & GRUPPIERUNG
# ==========================================
def get_plot_style(grp_name):
    grp_lower = grp_name.lower()
    
    if 'dns' in grp_lower or 'wt' in grp_lower: color = '#0072B2'
    elif 'mtrc' in grp_lower: color = '#E69F00'
    elif 'mtrab' in grp_lower or 'mtrb' in grp_lower: color = '#009E73'
    elif 'pdsa' in grp_lower: color = '#CC79A7'
    else: color = '#7f7f7f'
        
    if 'abts' in grp_lower: linestyle = '--'
    elif 'pms' in grp_lower: linestyle = ':'
    elif 'pes' in grp_lower: linestyle = '-.'
    elif 'fecn' in grp_lower: linestyle = (0, (3, 1, 1, 1))
    else: linestyle = '-'
        
    return color, linestyle

def plotte_gruppen(gruppen_liste, dateiname):
    aktive_gruppen = [g for g in gruppen_liste if g in gruppen_zuordnung]
    if not aktive_gruppen: return
        
    plt.figure(figsize=(10, 6))
    for g in aktive_gruppen:
        wells = gruppen_zuordnung[g]
        mittelwert = df[wells].mean(axis=1)
        std_abw = df[wells].std(axis=1)
        angezeigter_name = display_names.get(g, g)
        
        # Abruf von Farbe und Linienstil über den ursprünglichen Gruppennamen
        c, ls = get_plot_style(g)
        
        # Anwendung von Farbe (color) und Linienstil (linestyle)
        p = plt.plot(zeit, mittelwert, linewidth=2.5, color=c, linestyle=ls, label=f"{angezeigter_name} (n={len(wells)})")
        plt.fill_between(zeit, mittelwert - std_abw, mittelwert + std_abw, color=c, alpha=0.15)

    plt.xlabel('Time [h]', fontsize=12)
    plt.ylabel(r'Optical Density (OD$_{600}$)', fontsize=12)
    plt.xlim(0, 24)                     
    plt.xticks(range(0, 25, 4))
    plt.grid(True, linestyle=':', alpha=0.4)
    plt.legend(
    title="Mean ± SD",
    loc="lower center",
    bbox_to_anchor=(0.5, 1.02),
    ncol=3,
    frameon=False
)
    plt.tight_layout()
    plt.savefig(os.path.join(output_ordner, f"{datei_datum}{dateiname}"), dpi=360, bbox_inches='tight')
    plt.close()

# GRUPPEN FÜR DIE AUSWERTUNG DEFINIEREN
gruppen_ohne = [g for g in gruppen_zuordnung.keys() if g not in blank_gruppen and not any(m in g for m in mediatoren)]
gruppen_abts = [g for g in gruppen_zuordnung.keys() if 'ABTS' in g and g not in blank_gruppen]
gruppen_pms  = [g for g in gruppen_zuordnung.keys() if 'PMS' in g and g not in blank_gruppen]
gruppen_pes  = [g for g in gruppen_zuordnung.keys() if 'PES' in g and g not in blank_gruppen]
gruppen_fecn = [g for g in gruppen_zuordnung.keys() if 'FeCN' in g and g not in blank_gruppen]
alle_gruppen = gruppen_ohne + gruppen_abts + gruppen_pms + gruppen_pes + gruppen_fecn

plotte_gruppen(gruppen_ohne, 'Wachstum_ohne_Mediator.png')
plotte_gruppen(gruppen_abts, 'Wachstum_ABTS.png')
plotte_gruppen(gruppen_pms, 'Wachstum_PMS.png')
plotte_gruppen(gruppen_pes, 'Wachstum_PES.png')
plotte_gruppen(gruppen_fecn, 'Wachstum_FeCN.png')
plotte_gruppen(alle_gruppen, 'Wachstum_Alle_Mediatoren.png')

# ==========================================
# 6. TEXT-DATEI (KINETIK & T-TESTS MIT KORREKTUR)
# ==========================================
root_msg = tk.Tk()
root_msg.withdraw()
interaktiv_modus = messagebox.askyesno(
    "Interaktive Auswertung", 
    "Möchtest du die Wachstumsphasen für jede Gruppe grafisch überprüfen und interaktiv anpassen?\n\n(Ja = Ein Fenster pro Gruppe öffnet sich\nNein = Automatische Berechnung für alle)"
)
root_msg.destroy()

max_ods = {g: df[wells].max() for g, wells in gruppen_zuordnung.items()}
aucs = {g: df[wells].apply(lambda col: np.trapezoid(col.dropna(), zeit[col.notna()])) for g, wells in gruppen_zuordnung.items()}

def holm_bonferroni(p_values_dict, alpha=0.05):
    sorted_tests = sorted(p_values_dict.items(), key=lambda x: x[1])
    m = len(sorted_tests)
    results = {}
    for rank, (name, p_val) in enumerate(sorted_tests):
        adjusted_alpha = alpha / (m - rank)
        sig = "Signifikant *" if p_val < adjusted_alpha else "n.s."
        results[name] = (p_val, adjusted_alpha, sig)
    return results

pfad_statistik = os.path.join(output_ordner, f"{datei_datum}Statistik_T_Tests.txt")
with open(pfad_statistik, "w", encoding='utf-8') as f:
    
    # --- 1. KINETISCHE WACHSTUMSPARAMETER ---
    f.write("=== KINETISCHE WACHSTUMSPARAMETER ===\n")
    f.write(f"{'Gruppe':<35} | {'Max OD':<8} | {'Zeit (Max)':<12} | {'Lag-Phase':<12} | {'Doubling Time':<15} | {'AUC':<10}\n")
    f.write("-" * 105 + "\n")
    
    for g in alle_gruppen:
        if g not in gruppen_zuordnung: continue
        wells = gruppen_zuordnung[g]
        mittelwert = df[wells].mean(axis=1)
        angezeigter_name = display_names.get(g, g)
        
        if zeit.max() > 1.0:
            gueltige_zeit_mask = zeit >= 1.0
            max_idx = mittelwert[gueltige_zeit_mask].idxmax()
        else:
            max_idx = mittelwert.idxmax()
            
        max_od_val = mittelwert[max_idx]
        max_t = zeit[max_idx]
        
        time_str = f"{max_t:.2f} h"

        t_clean = zeit[mittelwert.notna()]
        od_clean = mittelwert[mittelwert.notna()]
        
        rec_lag, rec_td, rec_mu, rec_b, rec_t_window, rec_y0 = berechne_kinetik_empfehlung(t_clean, od_clean)
        
        if interaktiv_modus:
            fitter = InteraktiverFitter(t_clean.values, od_clean.values, angezeigter_name, rec_lag, rec_td, rec_mu, rec_b, rec_t_window, rec_y0)
            lag, td = fitter.run()
        else:
            lag, td = rec_lag, rec_td
            
        auc_mean = aucs[g].mean() if g in aucs else 0
        
        lag_str = f"{lag:.2f} h" if lag is not None else "N/A"
        td_str = f"{td:.2f} h" if td is not None else "N/A"
        
        f.write(f"{angezeigter_name:<35} | {max_od_val:<8.3f} | {time_str:<12} | {lag_str:<12} | {td_str:<15} | {auc_mean:<10.2f}\n")
            
    f.write("\n\n")

    # --- 2. STATISTIK: MUTANTE VS KONTROLLE (INNERHALB EINER BEDINGUNG) ---
    f.write("=== STATISTISCHE AUSWERTUNG I: MUTANTEN VS. KONTROLLE ===\n")
    f.write("Methode: Welch's T-Test mit anschließender Holm-Bonferroni Korrektur.\n\n")
    
    metriken = [("End-Yield (Max OD)", max_ods), ("Fläche unter der Kurve (AUC)", aucs)]
    
    test_bedingungen = [
        ("Bedingung: Ohne Mediator", gruppen_ohne, 'Vnat_dns'),
        ("Bedingung: + ABTS", gruppen_abts, 'Vnat_ABTS_dns'),
        ("Bedingung: + PMS", gruppen_pms, 'Vnat_PMS_dns'),
        ("Bedingung: + PES", gruppen_pes, 'Vnat_PES_dns'),
        ("Bedingung: + FeCN", gruppen_fecn, 'Vnat_FeCN_dns')
    ]
    
    for metrik_name, daten_dict in metriken:
        f.write(f"--- PARAMETER: {metrik_name} ---\n")
        
        for cond_titel, cond_gruppen, cond_control in test_bedingungen:
            if cond_control not in daten_dict: continue
                
            control_name = display_names.get(cond_control, cond_control)
            f.write(f"\n*) {cond_titel} (Vergleich vs. Kontrolle: {control_name})\n")
            
            raw_p_values = {}
            for g in cond_gruppen:
                if g != cond_control and g in daten_dict:
                    _, p_val = stats.ttest_ind(daten_dict[cond_control], daten_dict[g], equal_var=False, nan_policy='omit')
                    raw_p_values[g] = p_val
            
            if raw_p_values:
                korrigiert = holm_bonferroni(raw_p_values)
                for mutante, (p, adj_alpha, sig) in korrigiert.items():
                    mutante_name = display_names.get(mutante, mutante)
                    f.write(f"  {control_name} vs. {mutante_name:<30}: p-Wert = {p:.4f} (Holm-Limit: {adj_alpha:.4f}) -> {sig}\n")
            else:
                f.write("  Keine Mutantengruppen zum Vergleichen in dieser Bedingung gefunden.\n")
                
        f.write("\n" + "="*60 + "\n\n")

# --- 3. STATISTIK: KREUZVERGLEICH (OHNE VS MIT MEDIATOR PRO STAMM) ---
    f.write("=== STATISTISCHE AUSWERTUNG II: EINFLUSS DES MEDIATORS PRO STAMM ===\n")
    f.write("Methode: Welch's T-Test (unkorrigiert).\n")
    f.write("Hinweis: Da es sich hierbei um unabhängige, gezielte Einzelvergleiche desselben Stammes handelt, wurde auf eine multiple Testkorrektur verzichtet, um falsch-negative Ergebnisse (Typ-II-Fehler) zu vermeiden.\n\n")
    
    for metrik_name, daten_dict in metriken:
        f.write(f"--- PARAMETER: {metrik_name} ---\n")
        
        for med in mediatoren:
            f.write(f"\n*) Einfluss von {med} auf die jeweiligen Stämme\n")
            gefunden = False
            
            for base_grp in gruppen_ohne:
                # Rekonstruiere den Gruppennamen mit Mediator (z.B. Vnat_dns -> Vnat_FeCN_dns)
                if base_grp.startswith("Vnat_"):
                    parts = base_grp.split("_", 1)
                    if len(parts) == 2:
                        med_grp = f"Vnat_{med}_{parts[1]}"
                        
                        if base_grp in daten_dict and med_grp in daten_dict:
                            _, p_val = stats.ttest_ind(daten_dict[base_grp], daten_dict[med_grp], equal_var=False, nan_policy='omit')
                            
                            sig = "Signifikant *" if p_val < 0.05 else "n.s."
                            
                            base_name = display_names.get(base_grp, base_grp)
                            med_name = display_names.get(med_grp, med_grp)
                            
                            f.write(f"  {base_name} vs. {med_name:<30}: p-Wert = {p_val:.4f} -> {sig}\n")
                            gefunden = True
            
            if not gefunden:
                f.write("  Keine passenden Stammpaare für diesen Vergleich gefunden.\n")
                
        f.write("\n" + "="*60 + "\n\n")

print("Auswertung abgeschlossen! Diagramme und Statistik-Textdatei befinden sich im Ordner:", output_ordner)
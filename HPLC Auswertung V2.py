import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox
from typing import Dict, List, Optional, Tuple
import re
from datetime import datetime
import scipy.stats as stats
from itertools import combinations

DEFAULT_COMPOUNDS = ['Glucose', 'Pyruvate', 'Lactate', 'Formate', 'Acetate', 'Ethanol']
AMOUNT_SUFFIX = '|Amount'
RT_SUFFIX = '|RT'

class HPLCDataAnalyzer:
    """Hauptklasse zur Analyse von HPLC-Daten"""
    
    def __init__(self):
        self.data: Optional[pd.DataFrame] = None
        self.filepath: Optional[Path] = None
        self.compounds: List[str] = DEFAULT_COMPOUNDS.copy()
        self.selected_experiments: Dict[str, bool] = {}
        self.plot_settings: Dict = {
            'figsize': (10, 6),      
            'dpi': 600,              
            'style': 'default',      
            'font_size': 12,         
            'save_format': 'png'
        }
    
    def load_data(self, filepath: str) -> bool:
        """Lädt Excel-/CSV-Daten robust und erkennt die Header-Zeile automatisch."""
        try:
            self.filepath = Path(filepath)

            if not self.filepath.exists():
                print(f"✗ ERROR: Datei nicht gefunden: {self.filepath}")
                return False

            suffix = self.filepath.suffix.lower()

            try:
                if suffix == '.csv':
                    self.data = pd.read_csv(
                        self.filepath,
                        sep=None,
                        engine='python',
                        encoding='utf-8-sig'
                    )
                elif suffix == '.xlsx':
                    self.data = pd.read_excel(
                        self.filepath,
                        sheet_name=0,
                        engine='openpyxl',
                        header=None
                    )
                elif suffix == '.xls':
                    self.data = pd.read_excel(
                        self.filepath,
                        sheet_name=0,
                        engine='xlrd',
                        header=None
                    )
                else:
                    print(f"✗ ERROR: Nicht unterstütztes Dateiformat: {suffix}")
                    return False

            except PermissionError:
                print("✗ ERROR: Datei kann nicht gelesen werden. Bitte schließen Sie die Datei in Excel.")
                return False
            except ImportError:
                if suffix == '.xlsx':
                    print("✗ ERROR: Für .xlsx-Dateien fehlt 'openpyxl'.")
                    print("Installieren mit: python -m pip install openpyxl")
                elif suffix == '.xls':
                    print("✗ ERROR: Für .xls-Dateien fehlt 'xlrd'.")
                    print("Installieren mit: python -m pip install xlrd")
                return False
            except ValueError as ve:
                print(f"✗ ERROR: Datei konnte nicht als {suffix}-Datei gelesen werden: {ve}")
                return False

            if self.data is None or self.data.empty:
                print("✗ ERROR: Die Datei enthält keine Daten.")
                return False

            # HPLC-Exporte enthalten häufig mehrere Metadatenzeilen.
            # Daher wird die Zeile mit der Sample-Spalte automatisch gesucht.
            sample_header_row = None

            for row_idx in range(min(len(self.data), 20)):
                row_values = self.data.iloc[row_idx].astype(str).str.strip()
                if row_values.str.lower().eq('sample').any():
                    sample_header_row = row_idx
                    break

            if sample_header_row is not None:
                self.data.columns = [
                    str(col).strip() if pd.notna(col) else ''
                    for col in self.data.iloc[sample_header_row]
                ]
                self.data = self.data.iloc[sample_header_row + 1:].reset_index(drop=True)
            else:
                # Fallback: erste Zeile als Header verwenden.
                self.data.columns = [
                    str(col).strip() if pd.notna(col) else ''
                    for col in self.data.iloc[0]
                ]
                self.data = self.data.iloc[1:].reset_index(drop=True)

            raw_columns = [
                str(col).strip() if pd.notna(col) else ''
                for col in self.data.columns
            ]
            
            unique_cols = []
            seen = {}
            for col in raw_columns:
                if col in seen:
                    seen[col] += 1
                    # Bei Duplikaten eine Nummer anhängen (z.B. '' wird zu '_1')
                    unique_cols.append(f"{col}_{seen[col]}")
                else:
                    seen[col] = 0
                    unique_cols.append(col)
                    
            self.data.columns = unique_cols

            print(f"\n✓ Gefundene Spalten: {list(self.data.columns)}")

            if 'Sample' not in self.data.columns:
                sample_cols = [
                    col for col in self.data.columns
                    if 'sample' in str(col).lower()
                ]

                if sample_cols:
                    print(f"⚠ 'Sample' nicht gefunden, verwende: {sample_cols[0]}")
                    self.data.rename(
                        columns={sample_cols[0]: 'Sample'},
                        inplace=True
                    )
                else:
                    print("✗ ERROR: Keine Sample-Spalte gefunden!")
                    print(f"Verfügbare Spalten: {list(self.data.columns)}")
                    return False

            self.data = self.data.dropna(how='all').copy()
            self.data = self.data[self.data['Sample'].notna()].copy()
            self.data['Sample'] = (
                self.data['Sample'].astype(str).str.strip()
            )

            print(f"✓ Daten erfolgreich geladen: {len(self.data)} Zeilen")
            print(f"✓ Erste 5 Samples: {list(self.data['Sample'].head())}")

            return True

        except Exception as e:
            print(f"✗ Fehler beim Laden: {e}")
            import traceback
            traceback.print_exc()
            return False

    def parse_sample_name(self, sample_name: str) -> Dict:
        """Parsen des Sample-Namens in Komponenten"""
        result = {
            'original': sample_name,
            'experiment_type': None,
            'date': None,
            'time_h': None,    
            'medium': None,
            'mutant': None,
            'suffix': None,
            'is_biolector': False
        }
        
        sample_name = str(sample_name).strip()
        invalid_samples = ['nan', 'None', '', 'HeaderName', 'NumberOfRows', 
                          'NumberOfCol', 'NumberOfHead', 'Modified']
        if sample_name in invalid_samples:
            return result
        
        
        if '_Biolector_' in sample_name:
            result['is_biolector'] = True
            parts = sample_name.split('_')
            if len(parts) >= 2:
                result['date'] = f"{parts[0]}_{parts[1]}" if parts[0].isdigit() else parts[0]
            for medium in ['M1', 'FeCN', 'ABTS']:
                if medium in sample_name:
                    result['medium'] = medium
                    break
            if result['medium']:
                mut_start = sample_name.find(result['medium']) + len(result['medium']) + 1
                mut_part = sample_name[mut_start:]               
                mut_clean = re.sub(r'[\s_]*\d+$', '', mut_part).strip()
               
                result['mutant'] = mut_clean
            result['experiment_type'] = 'biolector'
                           
        elif sample_name.startswith('Ecuev'):
            result['experiment_type'] = 'ecuev'
            if '_dns' in sample_name: result['mutant'] = 'dns'
            elif '_pdsA' in sample_name: result['mutant'] = 'pdsA'
            elif '_mtrC' in sample_name: result['mutant'] = 'mtrC'
            elif '_mtrAB' in sample_name: result['mutant'] = 'mtrAB'
            
            if sample_name.endswith('_P'): result['suffix'] = 'P'
            elif sample_name.endswith('_RM'): result['suffix'] = 'RM'
            elif sample_name.endswith('_C'): result['suffix'] = 'C'
            else: result['suffix'] = 'main'
        
        elif re.match(r'^\d{6}_', sample_name):
            result['experiment_type'] = 'dated'
            parts = sample_name.split('_')
            
            if len(parts) >= 1:
                result['date'] = parts[0]
            
            # Zeit in Stunden suchen (endet auf 'h')
            for part in parts:
                if part.endswith('h') and part[:-1].isdigit():
                    result['time_h'] = int(part[:-1])
                    break
            
            for medium in ['M1', 'FeCN', 'ABTS']:
                if medium in sample_name:
                    result['medium'] = medium
                    break
        
        return result
    
    def add_compound(self, compound_name: str) -> bool:
        """Fügt neuen Stoff zur Analyse hinzu"""
        if compound_name not in self.compounds:
            self.compounds.append(compound_name)
            print(f"✓ Stoff hinzugefügt: {compound_name}")
            return True
        return False
    
    def remove_compound(self, compound_name: str) -> bool:
        """Entfernt Stoff aus der Analyse"""
        if compound_name in self.compounds:
            self.compounds.remove(compound_name)
            print(f"✓ Stoff entfernt: {compound_name}")
            return True
        return False
    
    def get_amount_column(self, compound: str) -> str:
        """Gibt den korrekten Spaltennamen für Amount zurück"""
        return f"{compound}{AMOUNT_SUFFIX}"
    
    def get_rt_column(self, compound: str) -> str:
        """Gibt den korrekten Spaltennamen für Retention Time zurück"""
        return f"{compound}{RT_SUFFIX}"
    
    def group_ecuev_data(self) -> Dict:
        """Gruppiert Ecuev-Experimente"""
        if self.data is None:
            return {}
        
        grouped = {
            'by_mutant': {},
            'controls': {}
        }
        
        ecuev_mask = self.data['Sample'].astype(str).str.startswith('Ecuev')
        ecuev_data = self.data[ecuev_mask].copy()
        
        print(f"\n📊 Ecuev-Proben gefunden: {len(ecuev_data)}")
        
        for idx, row in ecuev_data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            mutant = parsed['mutant']
            suffix = parsed['suffix']
            
            if mutant:
                if suffix in ['P', 'RM', 'C']:
                    control_key = f"{mutant}_{suffix}"
                    if control_key not in grouped['controls']:
                        grouped['controls'][control_key] = []
                    grouped['controls'][control_key].append(row)
                else:
                    if mutant not in grouped['by_mutant']:
                        grouped['by_mutant'][mutant] = []
                    grouped['by_mutant'][mutant].append(row)
        
        results = {
            'mutants_mean': {},
            'controls': {}
        }
        
        for mutant, rows in grouped['by_mutant'].items():
            df_group = pd.DataFrame(rows)
            mean_row = df_group.mean(numeric_only=True)
            results['mutants_mean'][mutant] = {
                'mean': mean_row,
                'n_samples': len(rows),
                'std': df_group.std(numeric_only=True),
                'raw_data': df_group 
            }
        
        for control, rows in grouped['controls'].items():
            df_group = pd.DataFrame(rows)
            results['controls'][control] = {
                'data': df_group,
                'n_samples': len(rows)
            }
        
        return results
    
    def group_dated_data(self) -> Dict:
        """Gruppiert datumsbasierte Experimente nach Zeit (h) und Medium"""
        if self.data is None:
            return {}
        
        grouped = {'by_time_medium': {}}
        
        for idx, row in self.data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            
            if parsed['experiment_type'] == 'dated':
                if parsed['time_h'] is None or not parsed['medium']:
                    continue
                key = f"{parsed['time_h']}_{parsed['medium']}"
                if key not in grouped['by_time_medium']:
                    grouped['by_time_medium'][key] = []
                grouped['by_time_medium'][key].append(row)
        
        results = {}
        for key, rows in grouped['by_time_medium'].items():
            df_group = pd.DataFrame(rows)
            mean_row = df_group.mean(numeric_only=True)
            
            time_str, medium = key.split('_')
            
            results[key] = {
                'mean': mean_row,
                'n_samples': len(rows),
                'std': df_group.std(numeric_only=True),
                'time_h': int(time_str),
                'medium': medium
            }
        
        return results

    def calculate_dated_yield_timecourse(self) -> Optional[pd.DataFrame]:
        """Berechnet kumulative Yields fuer jeden Zeitpunkt der datierten Experimente.

        Fuer jedes Replikat und jeden Zeitpunkt gilt:
            Y_p/s(t) = [c_Produkt(t) - c_Produkt(t0)] /
                       [c_Glucose(t0) - c_Glucose(t)]

        Damit wird fuer jeden Messpunkt genau die bis dahin verbrauchte Glucose
        verwendet. Der t=0-Punkt selbst besitzt definitionsgemaess keinen Yield
        (0/0) und wird deshalb als NaN ausgegeben.
        """
        if self.data is None or self.data.empty:
            return None

        glucose_col = self.get_amount_column('Glucose')
        if glucose_col not in self.data.columns:
            print("⚠ Konnte Glucose-Spalte fuer Yield-Zeitverlauf nicht finden.")
            return None

        records = []
        for idx, row in self.data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            if parsed['experiment_type'] != 'dated':
                continue
            if parsed['time_h'] is None or not parsed['medium']:
                continue

            # Replikat-ID robust aus dem letzten numerischen Namensbestandteil lesen,
            # z.B. 260707_24h_M1_2 -> Replicate 2.
            parts = str(row['Sample']).strip().split('_')
            replicate = None
            for part in reversed(parts):
                if part.isdigit():
                    replicate = int(part)
                    break
            if replicate is None:
                # Fallback: falls keine Replikatnummer vorhanden ist, Sample-Name nutzen.
                replicate = str(row['Sample'])

            rec = {
                'original_idx': idx,
                'Sample': row['Sample'],
                'Medium': parsed['medium'],
                'Time_h': parsed['time_h'],
                'Replicate': replicate,
            }
            for compound in self.compounds:
                col = self.get_amount_column(compound)
                if col in self.data.columns:
                    rec[compound] = pd.to_numeric(row.get(col, np.nan), errors='coerce')
            records.append(rec)

        if not records:
            return None

        raw = pd.DataFrame(records)
        result_rows = []

        for medium in raw['Medium'].dropna().unique():
            med = raw[raw['Medium'] == medium].copy()
            for replicate in med['Replicate'].unique():
                rep = med[med['Replicate'] == replicate].sort_values('Time_h')
                if rep.empty:
                    continue

                # Der frueheste vorhandene Zeitpunkt dieses Replikats ist t0.
                start = rep.iloc[0]
                start_time = start['Time_h']
                glc_start = start.get('Glucose', np.nan)
                if pd.isna(glc_start):
                    continue

                for _, current in rep.iterrows():
                    time_h = current['Time_h']
                    glc_current = current.get('Glucose', np.nan)
                    consumed_glc = glc_start - glc_current if pd.notna(glc_current) else np.nan

                    out = {
                        'Medium': medium,
                        'Time_h': time_h,
                        'Time_d': time_h / 24.0,
                        'Replicate': replicate,
                        'Start_Time_h': start_time,
                        'Glucose_Start_g_L': glc_start,
                        'Glucose_Current_g_L': glc_current,
                        'Glucose_Consumed_g_L': consumed_glc,
                    }

                    for compound in self.compounds:
                        if compound.lower() == 'glucose':
                            continue
                        if compound not in rep.columns:
                            continue

                        comp_start = start.get(compound, np.nan)
                        comp_current = current.get(compound, np.nan)
                        if pd.isna(comp_start) or pd.isna(comp_current):
                             yield_val = np.nan
                        elif time_h == start_time:
                            yield_val = 0.0
                        elif pd.isna(consumed_glc) or consumed_glc <= 0:
                             yield_val = np.nan
                        else:
                            delta_product = comp_current - comp_start
                            yield_val = delta_product / consumed_glc    

                        out[f'Yield_{compound}_g_g'] = yield_val

                    result_rows.append(out)

        if not result_rows:
            return None
        return pd.DataFrame(result_rows)

    def plot_dated_yield_timecourse(self, yield_time_df: pd.DataFrame,
                                    compound: str,
                                    output_dir: str = './plots') -> Optional[str]:
        """Plottet den Yield-Zeitverlauf optisch analog zum Konzentrations-Zeitverlauf."""
        if yield_time_df is None or yield_time_df.empty:
            return None

        yield_col = f'Yield_{compound}_g_g'
        if yield_col not in yield_time_df.columns:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        fig, ax = plt.subplots(figsize=self.plot_settings['figsize'])

        colors = {'M1': '#0072B2', 'FeCN': '#E69F00', 'ABTS': '#CC79A7'}
        labels = {'M1': 'M1', 'FeCN': '5 mM FeCN', 'ABTS': '100 μM ABTS'}

        for medium in ['M1', 'FeCN', 'ABTS']:
            med = yield_time_df[yield_time_df['Medium'] == medium].copy()
            if med.empty:
                continue

            # Nur definierte Yields plotten; t=0 ist normalerweise NaN (0/0).
            summary = (med.groupby('Time_h')[yield_col]
                         .agg(['mean', 'std', 'count'])
                         .reset_index()
                         .sort_values('Time_h'))
            summary = summary[summary['mean'].notna()]
            if summary.empty:
                continue

            x_vals_days = summary['Time_h'].to_numpy(dtype=float) / 24.0
            y_vals = summary['mean'].to_numpy(dtype=float)
            # Bei n=1 ist pandas-std NaN; fuer die Darstellung dann 0 verwenden.
            y_err = summary['std'].fillna(0).to_numpy(dtype=float)

            ax.errorbar(x_vals_days, y_vals, yerr=y_err,
                        capsize=4,
                        label=labels.get(medium, medium),
                        color=colors.get(medium, 'gray'),
                        marker='o', linestyle='-', linewidth=2)

        ax.set_xlabel('Time [d]', fontsize=self.plot_settings['font_size'])
        ax.set_ylabel(f'{compound} Yield $Y_{{p/s}}$ [g/g]',
                      fontsize=self.plot_settings['font_size'])

        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)

        plt.tight_layout()
        filepath = str(Path(output_dir) / f"{compound}_yield_timecourse.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'],
                    bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()

        print(f"✓ Yield-Zeitverlauf erstellt: {filepath}")
        return filepath

    def calculate_dated_yields(self) -> Optional[pd.DataFrame]:
        """Berechnet die Ausbeuten (Yields) für die datierten Experimente zwischen Start- und Endtag."""
        if self.data is None or self.data.empty:
            return None

        # 1. Alle datierten Proben identifizieren und Metadaten extrahieren
        parsed_list = []
        for idx, row in self.data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            if parsed['experiment_type'] == 'dated' and parsed['date'] and parsed['medium']:
                parsed_list.append({
                    'original_idx': idx,
                    'date': parsed['date'],
                    'medium': parsed['medium'],
                    'sample_name': row['Sample']
                })

        if not parsed_list:
            return None

        meta_df = pd.DataFrame(parsed_list)
        dates = sorted(meta_df['date'].unique())
        
        # Wir brauchen mindestens zwei Messpunkte (Tag 0 und Tag 7)
        if len(dates) < 2:
            return None 

        start_date = dates[0]   # Erster Tag (t=0)
        end_date = dates[-1]    # Letzter Tag (t=Ende)

        # Überprüfen, ob Glucose vorhanden ist (wird als Referenz benötigt)
        glucose_col = self.get_amount_column('Glucose')
        if glucose_col not in self.data.columns:
            print("⚠ Konnte Glucose-Spalte für Yield-Berechnung nicht finden.")
            return None

        yield_results = []

        # 2. Replikate pro Medium matchen und Yields berechnen
        for medium in meta_df['medium'].unique():
            # Hole die Indizes für Start und Ende, alphabetisch sortiert nach Sample-Namen (damit Replikate matchen)
            start_indices = meta_df[(meta_df['date'] == start_date) & (meta_df['medium'] == medium)].sort_values('sample_name')['original_idx'].tolist()
            end_indices = meta_df[(meta_df['date'] == end_date) & (meta_df['medium'] == medium)].sort_values('sample_name')['original_idx'].tolist()

            # Bestimme die Anzahl der Replikate (z.B. 3)
            n_reps = min(len(start_indices), len(end_indices))

            for rep_idx in range(n_reps):
                start_row = self.data.iloc[start_indices[rep_idx]]
                end_row = self.data.iloc[end_indices[rep_idx]]

                try:
                    glc_start = float(start_row.get(glucose_col, 0))
                    glc_end = float(end_row.get(glucose_col, 0))
                except (ValueError, TypeError):
                    continue

                # Delta Glucose = Verbrauch (Start - Ende)
                delta_glc = glc_start - glc_end

                if delta_glc <= 0:
                    continue  # Wenn keine Glucose verbraucht wurde, kann kein Yield berechnet werden

                row_dict = {
                    'Medium': medium,
                    'Replicate': rep_idx + 1,
                    'Start_Date': start_date,
                    'End_Date': end_date,
                    'Delta_Glucose_g_L': delta_glc
                }

                # Berechne Deltas und Yields für alle anderen Stoffe
                for compound in self.compounds:
                    if compound.lower() == 'glucose':
                        continue
                        
                    comp_col = self.get_amount_column(compound)
                    if comp_col in start_row.index and comp_col in end_row.index:
                        try:
                            comp_start = float(start_row.get(comp_col, 0))
                            comp_end = float(end_row.get(comp_col, 0))
                            
                            # Delta Produkt = Produktion (Ende - Start)
                            delta_comp = comp_end - comp_start
                            yield_val = delta_comp / delta_glc

                            row_dict[f'Delta_{compound}_g_L'] = delta_comp
                            row_dict[f'Yield_{compound}_g_g'] = yield_val
                        except (ValueError, TypeError):
                            pass

                yield_results.append(row_dict)

        if not yield_results:
            return None

        yield_df = pd.DataFrame(yield_results)

        # 3. Mean und Standardabweichung pro Medium anhängen
        summary_rows = []
        for medium in yield_df['Medium'].unique():
            med_df = yield_df[yield_df['Medium'] == medium]
            
            mean_dict = {'Medium': medium, 'Replicate': 'MEAN', 'Start_Date': '', 'End_Date': ''}
            std_dict = {'Medium': medium, 'Replicate': 'STD', 'Start_Date': '', 'End_Date': ''}

            for col in yield_df.columns:
                if col not in ['Medium', 'Replicate', 'Start_Date', 'End_Date']:
                    mean_dict[col] = med_df[col].mean()
                    std_dict[col] = med_df[col].std()

            summary_rows.extend([mean_dict, std_dict])

        yield_df = pd.concat([yield_df, pd.DataFrame(summary_rows)], ignore_index=True)
        return yield_df

    def calculate_biolector_yields(self, biolector_data: Dict) -> Optional[pd.DataFrame]:
        """Berechnet die Ausbeuten (Yields) für die Biolector-Experimente."""
        if not biolector_data:
            return None
            
        yield_results = []
        glc_col = self.get_amount_column('Glucose')
        start_glucose = 9.09
        
        for key, info in biolector_data.items():
            medium = info.get('medium')
            mutant = info.get('mutant')
            raw_df = info.get('raw_data')
            
            if raw_df is None or raw_df.empty or glc_col not in raw_df.columns:
                continue
                
            consumed_glc = start_glucose - raw_df[glc_col]
            
            row = {'Gruppe': key, 'Medium': medium, 'Mutante': mutant, 'n_Proben': len(raw_df)}
            
            for comp in self.compounds:
                if comp.lower() == 'glucose':
                    continue
                
                comp_col = self.get_amount_column(comp)
                if comp_col in raw_df.columns:
                    # Yields berechnen
                    yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                    valid_yields = yields[~np.isnan(yields)]
                    
                    if len(valid_yields) > 0:
                        mean_val = np.nanmean(valid_yields)
                        std_val = np.nanstd(valid_yields)
                        row[f'{comp}_Yield_mean'] = mean_val
                        row[f'{comp}_Yield_std'] = std_val
                        row[f'{comp} Yield (Mean ± SD)'] = f"{mean_val:.4f}±{std_val:.4f}"
                    else:
                        row[f'{comp}_Yield_mean'] = np.nan
                        row[f'{comp}_Yield_std'] = np.nan
                        row[f'{comp} Yield (Mean ± SD)'] = "N/A"
                        
            yield_results.append(row)
            
        if not yield_results:
            return None
            
        return pd.DataFrame(yield_results)

    def export_biolector_statistics(self, biolector_data: Dict, output_dir: str) -> Optional[str]:
        """Berechnet p-Werte (Welch's t-test) für die BioLector-Yields und exportiert sie als CSV."""
        if not biolector_data:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        stats_results = []
        glc_col = self.get_amount_column('Glucose')
        start_glucose = 9.09

        # Alle vorhandenen Medien und Mutanten auslesen
        medien = list(set([v.get('medium') for v in biolector_data.values() if v.get('medium')]))
        mutanten = list(set([v.get('mutant') for v in biolector_data.values() if v.get('mutant')]))
        
        plot_compounds = [c for c in self.compounds if c.lower() != 'glucose']

        # Hilfsfunktion zum sauberen Extrahieren der Yields
        def get_yields(medium, mutant, compound):
            key = f"{medium}_{mutant}"
            if key in biolector_data:
                raw_df = biolector_data[key]['raw_data']
                comp_col = self.get_amount_column(compound)
                if comp_col in raw_df.columns and glc_col in raw_df.columns:
                    consumed_glc = start_glucose - raw_df[glc_col]
                    yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                    return yields[~np.isnan(yields)]
            return np.array([])

        # 1. Mutanten-Vergleich (Referenz: dns) pro Medium
        if 'dns' in mutanten:
            for medium in medien:
                for comp in plot_compounds:
                    ctrl_data = get_yields(medium, 'dns', comp)
                    if len(ctrl_data) >= 2:
                        for mutant in mutanten:
                            if mutant == 'dns': continue
                            test_data = get_yields(medium, mutant, comp)
                            
                            if len(test_data) >= 2:
                                t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                                sig_level = '*** (p<0.001)' if p_val < 0.001 else '** (p<0.01)' if p_val < 0.01 else '* (p<0.05)' if p_val < 0.05 else 'ns'
                                
                                stats_results.append({
                                    'Vergleichsart': 'Mutanten (Ref: dns)',
                                    'Kategorie': f"Medium: {medium}",
                                    'Stoff': comp,
                                    'Vergleich': f"dns vs {mutant}",
                                    'n_Ref': len(ctrl_data),
                                    'n_Test': len(test_data),
                                    't_statistic': round(t_stat, 4),
                                    'p_value': round(p_val, 5),
                                    'Signifikanz': sig_level
                                })

        # 2. Medien-Vergleich (Referenz: M1) pro Mutante
        if 'M1' in medien:
            for mutant in mutanten:
                for comp in plot_compounds:
                    ctrl_data = get_yields('M1', mutant, comp)
                    if len(ctrl_data) >= 2:
                        for medium in medien:
                            if medium == 'M1': continue
                            test_data = get_yields(medium, mutant, comp)
                            
                            if len(test_data) >= 2:
                                t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                                sig_level = '*** (p<0.001)' if p_val < 0.001 else '** (p<0.01)' if p_val < 0.01 else '* (p<0.05)' if p_val < 0.05 else 'ns'
                                
                                stats_results.append({
                                    'Vergleichsart': 'Medien (Ref: M1)',
                                    'Kategorie': f"Mutante: {mutant}",
                                    'Stoff': comp,
                                    'Vergleich': f"M1 vs {medium}",
                                    'n_Ref': len(ctrl_data),
                                    'n_Test': len(test_data),
                                    't_statistic': round(t_stat, 4),
                                    'p_value': round(p_val, 5),
                                    'Signifikanz': sig_level
                                })

        if not stats_results:
            return None

        # Tabelle erstellen und exportieren
        stats_df = pd.DataFrame(stats_results)
        filepath = str(Path(output_dir) / "biolector_yield_statistics.csv")
        stats_df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
        return filepath

    def calculate_ecuev_yields(self, ecuev_data: Dict, start_glucose: float = 9.09) -> Optional[pd.DataFrame]:
        """Berechnet die Ausbeuten (Yields) für die E-Cuvette-Experimente."""
        if not ecuev_data or 'mutants_mean' not in ecuev_data:
            return None
            
        yield_results = []
        glc_col = self.get_amount_column('Glucose')
        
        # Gehe durch die Haupt-Mutanten
        for mutant, info in ecuev_data['mutants_mean'].items():
            raw_df = info.get('raw_data')
            
            if raw_df is None or raw_df.empty or glc_col not in raw_df.columns:
                continue
                
            consumed_glc = start_glucose - raw_df[glc_col]
            row = {'Gruppe': 'Ecuev', 'Mutante': mutant, 'n_Proben': len(raw_df)}
            
            for comp in self.compounds:
                if comp.lower() == 'glucose':
                    continue
                
                comp_col = self.get_amount_column(comp)
                if comp_col in raw_df.columns:
                    yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                    valid_yields = yields[~np.isnan(yields)]
                    
                    if len(valid_yields) > 0:
                        mean_val = np.nanmean(valid_yields)
                        std_val = np.nanstd(valid_yields)
                        row[f'{comp}_Yield_mean'] = mean_val
                        row[f'{comp}_Yield_std'] = std_val
                        row[f'{comp} Yield (Mean ± SD)'] = f"{mean_val:.4f}±{std_val:.4f}"
                    else:
                        row[f'{comp}_Yield_mean'] = np.nan
                        row[f'{comp}_Yield_std'] = np.nan
                        row[f'{comp} Yield (Mean ± SD)'] = "N/A"
                        
            yield_results.append(row)
            
        return pd.DataFrame(yield_results) if yield_results else None

    def export_ecuev_statistics(self, ecuev_data: Dict, output_dir: str, start_glucose: float = 9.09) -> Optional[str]:
        """Berechnet p-Werte (Welch-Test) für die Ecuev-Yields und exportiert sie als CSV."""
        if not ecuev_data or 'mutants_mean' not in ecuev_data:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        stats_results = []
        glc_col = self.get_amount_column('Glucose')
        plot_compounds = [c for c in self.compounds if c.lower() != 'glucose']
        mutants_data = ecuev_data['mutants_mean']

        def get_yields(mutant, compound):
            if mutant in mutants_data:
                raw_df = mutants_data[mutant].get('raw_data')
                comp_col = self.get_amount_column(compound)
                if raw_df is not None and comp_col in raw_df.columns and glc_col in raw_df.columns:
                    consumed_glc = start_glucose - raw_df[glc_col]
                    yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                    return yields[~np.isnan(yields)]
            return np.array([])

        mutanten = list(mutants_data.keys())
        
        # Mutanten-Vergleich (Referenz: dns)
        if 'dns' in mutanten:
            for comp in plot_compounds:
                ctrl_data = get_yields('dns', comp)
                if len(ctrl_data) >= 2:
                    for mutant in mutanten:
                        if mutant == 'dns': continue
                        test_data = get_yields(mutant, comp)
                        
                        if len(test_data) >= 2:
                            t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                            sig_level = '*** (p<0.001)' if p_val < 0.001 else '** (p<0.01)' if p_val < 0.01 else '* (p<0.05)' if p_val < 0.05 else 'ns'
                            
                            stats_results.append({
                                'Vergleichsart': 'E-Cuvette Mutanten (Ref: dns)',
                                'Stoff': comp,
                                'Vergleich': f"dns vs {mutant}",
                                'n_Ref': len(ctrl_data),
                                'n_Test': len(test_data),
                                't_statistic': round(t_stat, 4),
                                'p_value': round(p_val, 5),
                                'Signifikanz': sig_level
                            })

        if not stats_results:
            return None

        stats_df = pd.DataFrame(stats_results)
        filepath = str(Path(output_dir) / "ecuev_yield_statistics.csv")
        stats_df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
        return filepath

    def export_yield_statistics(self, yield_df: pd.DataFrame, output_dir: str) -> Optional[str]:
        """Berechnet p-Werte (Welch's t-test) für die Yield-Unterschiede zwischen den Medien."""
        if yield_df is None or yield_df.empty:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Filtere nur die echten Messwerte heraus (MEAN und STD Zeilen ignorieren)
        raw_df = yield_df[~yield_df['Replicate'].isin(['MEAN', 'STD'])].copy()
        
        media = raw_df['Medium'].unique()
        if len(media) < 2:
            return None  # Statistik braucht mindestens zwei Gruppen zum Vergleichen

        stats_results = []
        
        # Alle möglichen Paare bilden (z.B. M1 vs FeCN, M1 vs ABTS)
        media_pairs = list(combinations(media, 2))

        # Für jeden Stoff (außer Glucose) den Test durchführen
        for compound in self.compounds:
            if compound.lower() == 'glucose':
                continue
                
            yield_col = f'Yield_{compound}_g_g'
            
            if yield_col not in raw_df.columns:
                continue

            for m1, m2 in media_pairs:
                # Daten für die beiden Medien extrahieren
                data1 = raw_df[raw_df['Medium'] == m1][yield_col].dropna()
                data2 = raw_df[raw_df['Medium'] == m2][yield_col].dropna()

                # Wir brauchen mindestens 2 Werte pro Gruppe für einen t-Test
                if len(data1) >= 2 and len(data2) >= 2:
                    # Welch's t-test (equal_var=False)
                    t_stat, p_val = stats.ttest_ind(data1, data2, equal_var=False)
                    
                    # Signifikanz-Level bestimmen
                    if p_val < 0.001:
                        sig_level = '*** (p < 0.001)'
                    elif p_val < 0.01:
                        sig_level = '** (p < 0.01)'
                    elif p_val < 0.05:
                        sig_level = '* (p < 0.05)'
                    else:
                        sig_level = 'ns (not significant)'

                    stats_results.append({
                        'Stoff': compound,
                        'Vergleich': f"{m1} vs {m2}",
                        'Medium_1_n': len(data1),
                        'Medium_2_n': len(data2),
                        't_statistic': round(t_stat, 4),
                        'p_value': round(p_val, 5),
                        'Signifikanz': sig_level
                    })

        if not stats_results:
            return None

        # Ergebnisse als DataFrame speichern und exportieren
        stats_df = pd.DataFrame(stats_results)
        filepath = str(Path(output_dir) / "dated_yields_statistics.csv")
        stats_df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
        
        return filepath
    
    def export_dated_end_concentration_statistics(self, output_dir: str) -> Optional[str]:
        """
        Vergleicht die absoluten Endkonzentrationen der datierten/timed Experimente.

        Fuer jeden Stoff wird der letzte vorhandene Zeitpunkt (max. Time_h) verwendet.
        M1 dient als Referenz und wird separat mit FeCN bzw. ABTS verglichen.
        Statistik: zweiseitiger ungepaarter Welch-t-Test (equal_var=False).
        """
        if self.data is None or self.data.empty:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        # Alle datierten Proben inklusive Zeit und Medium sammeln.
        records = []
        for _, row in self.data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            if parsed['experiment_type'] != 'dated':
                continue
            if parsed['time_h'] is None or not parsed['medium']:
                continue

            rec = {
                'Sample': row['Sample'],
                'Time_h': parsed['time_h'],
                'Medium': parsed['medium'],
            }

            for compound in self.compounds:
                col = self.get_amount_column(compound)
                if col in self.data.columns:
                    rec[compound] = pd.to_numeric(row.get(col, np.nan), errors='coerce')

            records.append(rec)

        if not records:
            return None

        dated_raw = pd.DataFrame(records)
        end_time_h = dated_raw['Time_h'].max()
        end_df = dated_raw[dated_raw['Time_h'] == end_time_h].copy()

        if end_df.empty or 'M1' not in end_df['Medium'].unique():
            return None

        test_media = [m for m in ['FeCN', 'ABTS'] if m in end_df['Medium'].unique()]
        if not test_media:
            return None

        stats_results = []

        # Auch Glucose einschliessen, da es sich hier um Endkonzentrationen handelt.
        for compound in sorted(self.compounds, key=str.lower):
            if compound not in end_df.columns:
                continue

            control_data = pd.to_numeric(
                end_df[end_df['Medium'] == 'M1'][compound], errors='coerce'
            ).dropna()

            if len(control_data) < 2:
                continue

            m1_mean = control_data.mean()
            m1_std = control_data.std(ddof=1)

            for medium in test_media:
                test_data = pd.to_numeric(
                    end_df[end_df['Medium'] == medium][compound], errors='coerce'
                ).dropna()

                if len(test_data) < 2:
                    continue

                test_mean = test_data.mean()
                test_std = test_data.std(ddof=1)
                difference = test_mean - m1_mean

                if pd.notna(m1_mean) and m1_mean != 0:
                    difference_percent = (difference / m1_mean) * 100.0
                else:
                    difference_percent = np.nan

                t_stat, p_val = stats.ttest_ind(
                    control_data, test_data, equal_var=False, nan_policy='omit'
                )

                if p_val < 0.001:
                    sig_level = '*** (p < 0.001)'
                elif p_val < 0.01:
                    sig_level = '** (p < 0.01)'
                elif p_val < 0.05:
                    sig_level = '* (p < 0.05)'
                else:
                    sig_level = 'ns (not significant)'

                stats_results.append({
                    'End_Time_h': end_time_h,
                    'Stoff': compound,
                    'Vergleich': f'M1 vs {medium}',
                    'M1_n': len(control_data),
                    'M1_Mean_g_L': round(m1_mean, 4),
                    'M1_SD_g_L': round(m1_std, 4),
                    f'{medium}_n': len(test_data),
                    f'{medium}_Mean_g_L': round(test_mean, 4),
                    f'{medium}_SD_g_L': round(test_std, 4),
                    'Difference_Test_minus_M1_g_L': round(difference, 4),
                    'Difference_percent_vs_M1': round(difference_percent, 2) if pd.notna(difference_percent) else np.nan,
                    't_statistic': round(t_stat, 4),
                    'p_value': round(p_val, 5),
                    'Signifikanz': sig_level,
                })

        if not stats_results:
            return None

        stats_df = pd.DataFrame(stats_results)
        filepath = str(Path(output_dir) / 'dated_end_concentration_statistics.csv')
        stats_df.to_csv(
            filepath,
            index=False,
            sep=';',
            decimal=',',
            encoding='utf-8-sig'
        )
        return filepath

    def group_biolector_data(self) -> Dict:
        """Gruppiert Biolector-Experimente nach Medium und Mutante"""
        if self.data is None:
            return {}
        
        grouped = {'by_medium_mutant': {}}
        
        for idx, row in self.data.iterrows():
            parsed = self.parse_sample_name(row['Sample'])
            
            if parsed['is_biolector']:
                if not parsed['medium'] or not parsed['mutant']:
                    continue
                key = f"{parsed['medium']}_{parsed['mutant']}"
                if key not in grouped['by_medium_mutant']:
                    grouped['by_medium_mutant'][key] = []
                grouped['by_medium_mutant'][key].append(row)
        
        results = {}
        for key, rows in grouped['by_medium_mutant'].items():
            df_group = pd.DataFrame(rows)
            mean_row = df_group.mean(numeric_only=True)
            results[key] = {
                'mean': mean_row,
                'n_samples': len(rows),
                'std': df_group.std(numeric_only=True),
                'medium': key.split('_')[0] if '_' in key else None,
                'mutant': key.split('_')[1] if '_' in key else None,
                'raw_data': df_group  
            }
        
        return results
    
    def create_summary_table(self, data_dict: Dict, title: str) -> pd.DataFrame:
        """Erstellt Übersichtstabelle"""
        rows = []
        
        for key, info in data_dict.items():
            row = {'Gruppe': key, 'n_Proben': info.get('n_samples', 0)}
            
            for compound in self.compounds:
                amount_col = self.get_amount_column(compound)
                if amount_col in info.get('mean', pd.Series()).index:
                    mean_val = info['mean'][amount_col]
                    std_val = info.get('std', pd.Series()).get(amount_col, 0)
                    row[f'{compound} (Mean ± SD)'] = f"{mean_val:.4f}±{std_val:.4f}"
                    row[f'{compound}_mean'] = mean_val
                    row[f'{compound}_std'] = std_val
            
            rows.append(row)
        
        return pd.DataFrame(rows)
    
    def plot_compound_timecourse(self, data_dict: Dict, compound: str, 
                                 output_dir: str = './plots') -> Optional[str]:
        """Plottet Zeitverlauf eines Stoffes (X-Achse in Tagen)"""
        if not data_dict:
            return None
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        fig, ax = plt.subplots(figsize=self.plot_settings['figsize'])
        
        colors = {'M1': '#0072B2', 'FeCN': '#E69F00', 'ABTS': '#CC79A7'}
        labels = {'M1': 'M1', 'FeCN': '5 mM FeCN', 'ABTS': '100 \u03BCM ABTS'}
        
        # Alle vorhandenen Zeiten (in Stunden) extrahieren und aufsteigend sortieren
        times_h = sorted(list(set([info['time_h'] for info in data_dict.values() if 'time_h' in info])))
        
        for medium in ['M1', 'FeCN', 'ABTS']:
            x_vals_h = []
            y_vals = []
            y_err = []
            
            for time_h in times_h:
                key = f"{time_h}_{medium}"
                if key in data_dict:
                    info = data_dict[key]
                    amount_col = self.get_amount_column(compound)
                    if amount_col in info['mean'].index:
                        x_vals_h.append(time_h)
                        y_vals.append(info['mean'][amount_col])
                        y_err.append(info.get('std', pd.Series()).get(amount_col, 0))
            
            if x_vals_h:
                # Stunden in Tage umrechnen für korrekte proportionale Abstände auf der x-Achse
                x_vals_days = [h / 24.0 for h in x_vals_h]
                
                ax.errorbar(x_vals_days, y_vals, yerr=y_err, 
                           capsize=4,                                
                           label=labels.get(medium, medium),         
                           color=colors.get(medium, 'gray'),         
                           marker='o', linestyle='-', linewidth=2)
        
        ax.set_xlabel('Time [d]', fontsize=self.plot_settings['font_size'])
        ax.set_ylabel(f'{compound} [g/L]', fontsize=self.plot_settings['font_size'])
        
        
        # X-Achse so einstellen, dass bevorzugt ganze Tage oder halbe Tage (je nach Länge) angezeigt werden
        from matplotlib.ticker import MaxNLocator
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        
        ax.legend(loc='upper right')
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        filepath = str(Path(output_dir) / f"{compound}_timecourse.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'], 
                   bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        
        print(f"✓ Plot erstellt: {filepath}")
        return filepath

    def plot_all_compounds_by_medium(self, data_dict: Dict, medium: str,
                                     output_dir: str = './plots') -> Optional[str]:
        """
        Plottet den zeitlichen Verlauf aller Compounds für EIN Medium.
        Glucose wird auf der rechten y-Achse dargestellt,
        alle übrigen Compounds auf der linken y-Achse.
        """
        if not data_dict:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        fig, ax_left = plt.subplots(figsize=self.plot_settings['figsize'])
        ax_right = ax_left.twinx()

        # Okabe-Ito-Palette
        compound_colors = {
            'Glucose':  '#000000',  # Schwarz
            'Pyruvate': '#E69F00',  # Orange
            'Lactate':  '#56B4E9',  # Hellblau
            'Formate':  '#CC79A7',  # Grün
            'Acetate':  '#009E73',  # Magenta
            'Ethanol':  '#D55E00',  # Vermillion / Rot
        }

        fallback_colors = ['#D55E00', '#F0E442', '#999999']

        times_h = sorted(list(set(
            info['time_h']
            for info in data_dict.values()
            if 'time_h' in info
        )))

        fallback_index = 0
        lines_left = []
        lines_right = []

        # Produkte alphabetisch, Glucose bewusst am Ende
        ordered_compounds = sorted(
            [c for c in self.compounds if c.lower() != 'glucose'],
            key=str.lower
        )
        glucose_compounds = [c for c in self.compounds if c.lower() == 'glucose']

        for compound in ordered_compounds + glucose_compounds:
            amount_col = self.get_amount_column(compound)

            x_vals_h = []
            y_vals = []
            y_err = []

            for time_h in times_h:
                key = f"{time_h}_{medium}"
                if key not in data_dict:
                    continue

                info = data_dict[key]

                if amount_col not in info['mean'].index:
                    continue

                mean_val = info['mean'][amount_col]
                std_val = info.get('std', pd.Series()).get(amount_col, 0)

                if pd.isna(mean_val):
                    continue

                x_vals_h.append(time_h)
                y_vals.append(mean_val)
                y_err.append(0 if pd.isna(std_val) else std_val)

            if not x_vals_h:
                continue

            x_vals_days = [h / 24.0 for h in x_vals_h]

            color = compound_colors.get(compound)
            if color is None:
                color = fallback_colors[fallback_index % len(fallback_colors)]
                fallback_index += 1

            # Glucose auf rechte Achse, alles andere links
            target_ax = ax_right if compound.lower() == 'glucose' else ax_left

            container = target_ax.errorbar(
                x_vals_days,
                y_vals,
                yerr=y_err,
                capsize=4,
                label=compound,
                color=color,
                marker='o',
                linestyle='-',
                linewidth=2
            )

            if compound.lower() == 'glucose':
                lines_right.append(container)
            else:
                lines_left.append(container)

        ax_left.set_xlabel('Time [d]', fontsize=self.plot_settings['font_size'])
        ax_left.set_ylabel('Products [g/L]', fontsize=self.plot_settings['font_size'])
        ax_right.set_ylabel('Glucose [g/L]', fontsize=self.plot_settings['font_size'])

        from matplotlib.ticker import MaxNLocator
        ax_left.xaxis.set_major_locator(MaxNLocator(integer=True))

        # Deine gewünschten Achsenbereiche:
        ax_left.set_ylim(0, 3)
        ax_right.set_ylim(0, 10)

        ax_left.grid(True, alpha=0.3)

        # Gemeinsame Legende
        handles = lines_left + lines_right
        labels = [h.get_label() for h in handles]
        ax_left.legend(handles, labels, loc='upper right')

        plt.tight_layout()

        filepath = str(Path(output_dir) / f"{medium}_all_compounds_timecourse.png")

        plt.savefig(
            filepath,
            dpi=self.plot_settings['dpi'],
            bbox_inches='tight',
            format=self.plot_settings['save_format']
        )
        plt.close()

        print(f"✓ Alle Compounds für {medium} geplottet: {filepath}")
        return filepath   

    
    def plot_medium_comparison(self, data_dict: Dict, compound: str,
                               output_dir: str = './plots') -> Optional[str]:
        """Vergleicht alle 3 Medien für einen Stoff"""
        if not data_dict:
            return None
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        medium_data = {'M1': [], 'FeCN': [], 'ABTS': []}
        
        for key, info in data_dict.items():
            medium = info.get('medium')
            if medium in medium_data:
                amount_col = self.get_amount_column(compound)
                if amount_col in info['mean'].index:
                    medium_data[medium].append(info['mean'][amount_col])
        
        fig, ax = plt.subplots(figsize=(10, 6))
        
        keys = ['M1', 'FeCN', 'ABTS']
        display_labels = ['M1', '5 mM FeCN', '100 \u03BCM ABTS']
        
        values = [np.mean(medium_data[m]) if medium_data[m] else 0 for m in keys]
        stds = [np.std(medium_data[m]) if medium_data[m] else 0 for m in keys]
        
        # Geplottet werden jetzt die ausführlichen Namen (display_labels)
        bars = ax.bar(display_labels, values, yerr=stds, capsize=5,
                     color=['#2ecc71', '#e74c3c', '#3498db'],
                     edgecolor='black', linewidth=1.5, alpha=0.8)
        
        ax.set_ylabel(f'{compound} [g/L] (Mean ± SD)', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        
        for bar, val, std in zip(bars, values, stds):
            height = bar.get_height()
            ax.annotate(f'{val:.3f}\n±{std:.3f}',
                       xy=(bar.get_x() + bar.get_width() / 2, height + std),
                       ha='center', va='bottom', 
                       fontsize=self.plot_settings['font_size']-1)
        
        plt.tight_layout()
        filepath = str(Path(output_dir) / f"{compound}_medium_comparison.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'],
                   bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        
        print(f"✓ Plot erstellt: {filepath}")
        return filepath

    def plot_dated_yields(self, yield_df: pd.DataFrame, output_dir: str = './plots') -> List[str]:
        """Plottet die berechneten Yields als Balkendiagramm mit Fehlerbalken und Signifikanzsternen."""
        if yield_df is None or yield_df.empty:
            return []

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        generated_plots = []

        # Filtere die vorberechneten Mean- und Std-Werte aus der Yield-Tabelle
        mean_df = yield_df[yield_df['Replicate'] == 'MEAN'].set_index('Medium')
        std_df = yield_df[yield_df['Replicate'] == 'STD'].set_index('Medium')
        
        # Rohdaten für den Signifikanztest extrahieren
        raw_df = yield_df[~yield_df['Replicate'].isin(['MEAN', 'STD'])].copy()

        # Definiere die Reihenfolge und Farben der Medien
        media_order = ['M1', 'FeCN', 'ABTS']
        plot_media = [m for m in media_order if m in mean_df.index]
        if not plot_media:
            plot_media = list(mean_df.index)

        # Zuordnung für die professionellen Achsen-Beschriftungen
        label_map = {'M1': 'M1', 'FeCN': '5 mM FeCN', 'ABTS': '100 \u03BCM ABTS'}
        display_labels = [label_map.get(m, m) for m in plot_media]

        colors = {'M1': '#0072B2', 'FeCN': '#E69F00', 'ABTS': '#CC79A7'}
        bar_colors = [colors.get(m, 'gray') for m in plot_media]

        # Erstelle für jeden Stoff (außer Glucose) einen eigenen Yield-Plot in alphabetischer Reihenfolge
        for compound in sorted(self.compounds):
            if compound.lower() == 'glucose':
                continue

            yield_col = f'Yield_{compound}_g_g'
            
            # Prüfen, ob der Stoff erfolgreich berechnet wurde
            if yield_col in mean_df.columns and not mean_df[yield_col].isna().all():
                fig, ax = plt.subplots(figsize=self.plot_settings['figsize'])
                
                values = [mean_df.loc[m, yield_col] for m in plot_media]
                stds = [std_df.loc[m, yield_col] for m in plot_media]
                
                # Hier werden nun die display_labels für die x-Achse übergeben
                bars = ax.bar(display_labels, values, yerr=stds, capsize=5,
                             color=bar_colors, edgecolor='black', linewidth=1.5, alpha=0.8)
                
                ax.set_ylabel(f'{compound} Yield $Y_{{p/s}}$ [g/g]', 
                            fontsize=self.plot_settings['font_size'])
                ax.grid(True, alpha=0.3, axis='y')
                
                # Maximale Höhe für die Signifikanzklammern ermitteln
                max_bar_height = 0
                
                # Werte und Standardabweichung über den Balken anzeigen
                for bar, val, std in zip(bars, values, stds):
                    height = bar.get_height()
                    if not pd.isna(val) and not pd.isna(std):
                        y_pos = max(0, height) + max(0, std)
                        if y_pos > max_bar_height:
                            max_bar_height = y_pos
                            
                        ax.annotate(f'{val:.3f}\n±{std:.3f}',
                                   xy=(bar.get_x() + bar.get_width() / 2, y_pos),
                                   ha='center', va='bottom', 
                                   fontsize=self.plot_settings['font_size']-1)
                
                control_medium = 'M1'
                if control_medium in plot_media:
                    ctrl_idx = plot_media.index(control_medium)
                    control_data = raw_df[raw_df['Medium'] == control_medium][yield_col].dropna()
                    
                    if len(control_data) >= 2:
                        # Abstandshalter für die Klammern
                        y_offset = max_bar_height * 0.15  
                        current_h = max_bar_height + y_offset
                        
                        for i, medium in enumerate(plot_media):
                            if medium == control_medium:
                                continue
                                
                            test_data = raw_df[raw_df['Medium'] == medium][yield_col].dropna()
                            if len(test_data) >= 2:
                                # Welch's t-test berechnen: zweistichprobig ungepaart zweiseitig für ungleiche varianzen
                                t_stat, p_val = stats.ttest_ind(control_data, test_data, equal_var=False)
                                
                                # Sterne zuweisen
                                if p_val < 0.05:
                                    if p_val < 0.001: text = '***'
                                    elif p_val < 0.01: text = '**'
                                    else: text = '*'
                                    
                                    # Koordinaten für die Klammer
                                    x1, x2 = ctrl_idx, i
                                    # Klammer-Höhe
                                    h_tick = max_bar_height * 0.02
                                    
                                    # Linien zeichnen
                                    ax.plot([x1, x1, x2, x2], [current_h, current_h+h_tick, current_h+h_tick, current_h], lw=1.2, c='black')
                                    # Sternchen setzen
                                    ax.text((x1+x2)*.5, current_h+h_tick, text, ha='center', va='bottom', color='black', fontsize=self.plot_settings['font_size'] + 2)
                                    
                                    # y-Ebene für die nächste mögliche Klammer anheben, damit sie nicht überlappen
                                    current_h += (max_bar_height * 0.12)
                        
                        # Y-Achse nach oben anpassen, damit die Klammern nicht abgeschnitten werden
                        ax.set_ylim(0, current_h + (max_bar_height * 0.05))
                # =========================================================

                plt.tight_layout()
                filepath = str(Path(output_dir) / f"{compound}_yield_comparison.png")
                plt.savefig(filepath, dpi=self.plot_settings['dpi'],
                           bbox_inches='tight', format=self.plot_settings['save_format'])
                plt.close()
                
                generated_plots.append(filepath)
                print(f"✓ Yield-Plot mit Signifikanzen erstellt: {filepath}")

        return generated_plots

    def plot_dated_yields_m1_fecn_all_compounds(self, yield_df: pd.DataFrame,
                                                 output_dir: str = './plots') -> Optional[str]:
        """
        Gemeinsamer Yield-Balkenplot fuer alle Produkte der Timed Data.
        Pro Compound stehen M1 und FeCN direkt nebeneinander; die Compound-Farbe
        bleibt identisch, FeCN wird zusaetzlich durch Schraffur gekennzeichnet.
        Signifikanz: zweiseitiger ungepaarter Welch-t-Test M1 vs FeCN je Compound.
        """
        if yield_df is None or yield_df.empty:
            return None

        Path(output_dir).mkdir(parents=True, exist_ok=True)

        raw_df = yield_df[~yield_df['Replicate'].isin(['MEAN', 'STD'])].copy()
        mean_df = yield_df[yield_df['Replicate'] == 'MEAN'].set_index('Medium')
        std_df = yield_df[yield_df['Replicate'] == 'STD'].set_index('Medium')

        if 'M1' not in mean_df.index or 'FeCN' not in mean_df.index:
            return None

        compound_colors = {
            'Pyruvate': '#E69F00',
            'Lactate':  '#56B4E9',
            'Formate':  '#CC79A7',
            'Acetate':  '#009E73',
            'Ethanol':  '#D55E00',
        }

        compounds = sorted(
            [c for c in self.compounds if c.lower() != 'glucose'
             and f'Yield_{c}_g_g' in mean_df.columns],
            key=str.lower
        )
        if not compounds:
            return None

        fig, ax = plt.subplots(figsize=(12, 6))
        x = np.arange(len(compounds), dtype=float)
        width = 0.34

        m1_vals, fecn_vals = [], []
        m1_stds, fecn_stds = [], []
        for compound in compounds:
            col = f'Yield_{compound}_g_g'
            m1_vals.append(mean_df.loc['M1', col])
            fecn_vals.append(mean_df.loc['FeCN', col])
            m1_stds.append(std_df.loc['M1', col])
            fecn_stds.append(std_df.loc['FeCN', col])

        # Gleiche Compound-Farbe fuer beide Medien; FeCN wird durch Schraffur unterschieden.
        colors = [compound_colors.get(c, '#999999') for c in compounds]
        bars_m1 = ax.bar(
            x - width / 2, m1_vals, width, yerr=m1_stds, capsize=5,
            color=colors, edgecolor='black', linewidth=1.2, alpha=0.8, label='M1'
        )
        bars_fecn = ax.bar(
            x + width / 2, fecn_vals, width, yerr=fecn_stds, capsize=5,
            color=colors, edgecolor='black', linewidth=1.2, alpha=0.8,
            hatch='//', label='5 mM FeCN'
        )

        ax.set_xticks(x)
        ax.set_xticklabels(compounds)
        ax.set_ylabel(r'Yield $Y_{p/s}$ [g/g]', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        ax.legend(loc='upper right')

        # Hoechsten sichtbaren Balken inkl. SD als Basis fuer Signifikanzklammern nutzen.
        tops = []
        for val, sd in zip(m1_vals + fecn_vals, m1_stds + fecn_stds):
            if pd.notna(val):
                tops.append(max(0, val) + (0 if pd.isna(sd) else max(0, sd)))
        global_max = max(tops) if tops else 1.0
        if global_max <= 0:
            global_max = 1.0

        highest_annotation = global_max
        for i, compound in enumerate(compounds):
            col = f'Yield_{compound}_g_g'
            data_m1 = raw_df[raw_df['Medium'] == 'M1'][col].dropna()
            data_fecn = raw_df[raw_df['Medium'] == 'FeCN'][col].dropna()

            if len(data_m1) < 2 or len(data_fecn) < 2:
                continue

            _, p_val = stats.ttest_ind(data_m1, data_fecn, equal_var=False)
            if p_val >= 0.05:
                continue

            stars = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'

            local_top = max(
                (0 if pd.isna(m1_vals[i]) else m1_vals[i]) +
                (0 if pd.isna(m1_stds[i]) else m1_stds[i]),
                (0 if pd.isna(fecn_vals[i]) else fecn_vals[i]) +
                (0 if pd.isna(fecn_stds[i]) else fecn_stds[i])
            )
            y = local_top + global_max * 0.08
            tick = global_max * 0.025
            x1, x2 = x[i] - width / 2, x[i] + width / 2
            ax.plot([x1, x1, x2, x2], [y, y + tick, y + tick, y],
                    lw=1.2, c='black')
            ax.text((x1 + x2) / 2, y + tick, stars, ha='center', va='bottom',
                    fontsize=self.plot_settings['font_size'] + 2)
            highest_annotation = max(highest_annotation, y + tick + global_max * 0.06)

        ax.set_ylim(bottom=0, top=highest_annotation + global_max * 0.08)
        plt.tight_layout()

        filepath = str(Path(output_dir) / 'all_compounds_yield_M1_vs_FeCN.png')
        plt.savefig(filepath, dpi=self.plot_settings['dpi'],
                    bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()

        print(f"✓ Gemeinsamer M1-vs-FeCN-Yield-Plot erstellt: {filepath}")
        return filepath
    
    def plot_biolector_mutant_comparison(self, biolector_data: Dict, 
                                          compound: str, medium: str,
                                          output_dir: str = './plots') -> Optional[str]:
        """Vergleicht Mutanten innerhalb eines Mediums für Biolector inkl. Signifikanz (Ref: dns)"""
        if not biolector_data:
            return None
        
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        medium_mutants = {k: v for k, v in biolector_data.items() if v.get('medium') == medium}
        if not medium_mutants:
            return None
            
        fig, ax = plt.subplots(figsize=self.plot_settings['figsize'])
        
        #Sortierung (dns -> pdsA -> mtrC -> mtrAB)
        custom_order = ['dns', 'pdsA', 'mtrC', 'mtrAB']
        mutants_raw = list(set(v['mutant'] for v in medium_mutants.values()))
        mutants = sorted(mutants_raw, key=lambda x: custom_order.index(x) if x in custom_order else 99)
            
        values = []
        stds = []
        amount_col = self.get_amount_column(compound)
        
        for m in mutants:
            key = f"{medium}_{m}"
            if key in medium_mutants and amount_col in medium_mutants[key]['mean'].index:
                values.append(medium_mutants[key]['mean'][amount_col])
                stds.append(medium_mutants[key].get('std', pd.Series()).get(amount_col, 0))
            else:
                values.append(0)
                stds.append(0)
                
        x_pos = np.arange(len(mutants))
        bars = ax.bar(x_pos, values, yerr=stds, capsize=5,
                     color=plt.cm.viridis(np.linspace(0.2, 0.8, len(mutants))),
                     edgecolor='black', linewidth=1.2, alpha=0.85)
        
        ax.set_xticks(x_pos)
        
        display_labels = [r'$\Delta\mathit{' + m + '}$' for m in mutants]
        ax.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=self.plot_settings['font_size'] + 2)
        
        ax.set_ylabel(f'{compound} [g/L] (Mean ± SD)', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        
        max_bar_height = max([v + s for v, s in zip(values, stds) if pd.notna(v) and pd.notna(s)] + [0])
        
        # Signifikanztest gegen 'dns'
        if 'dns' in mutants and amount_col:
            ctrl_key = f"{medium}_dns"
            if ctrl_key in medium_mutants:
                ctrl_data = medium_mutants[ctrl_key]['raw_data'][amount_col].dropna()
                ctrl_idx = mutants.index('dns')
                
                if len(ctrl_data) >= 2:
                    y_offset = max_bar_height * 0.15
                    current_h = max_bar_height + y_offset
                    
                    for i, m in enumerate(mutants):
                        if m == 'dns': continue
                        test_key = f"{medium}_{m}"
                        if test_key in medium_mutants:
                            test_data = medium_mutants[test_key]['raw_data'][amount_col].dropna()
                            if len(test_data) >= 2:
                                t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                                if p_val < 0.05:
                                    text = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
                                    x1, x2 = ctrl_idx, i
                                    h_tick = max_bar_height * 0.02
                                    ax.plot([x1, x1, x2, x2], [current_h, current_h+h_tick, current_h+h_tick, current_h], lw=1.2, c='black')
                                    ax.text((x1+x2)*.5, current_h+h_tick, text, ha='center', va='bottom', color='black', fontsize=self.plot_settings['font_size'] + 2)
                                    current_h += (max_bar_height * 0.12)
                    ax.set_ylim(0, current_h + (max_bar_height * 0.05))

        plt.tight_layout()
        filepath = str(Path(output_dir) / f"biolector_{medium}_{compound}_mutants.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'], bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        print(f"✓ Plot erstellt: {filepath}")
        return filepath

    def plot_biolector_yield_all_compounds(self, biolector_data: Dict, medium: str, output_dir: str = './plots') -> Optional[str]:
        """Plottet alle Compounds (Yield) für alle Mutanten in einem bestimmten Medium inkl. Signifikanz (Ref: dns)."""
        if not biolector_data:
            return None
            
        medium_mutants = {k: v for k, v in biolector_data.items() if v.get('medium') == medium}
        if not medium_mutants:
            return None
            
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # Sortierung
        custom_order = ['dns', 'pdsA', 'mtrC', 'mtrAB']
        mutants_raw = list(set(v['mutant'] for v in medium_mutants.values()))
        mutants = sorted(mutants_raw, key=lambda x: custom_order.index(x) if x in custom_order else 99)
            
        # ALPHABETISCHE SORTIERUNG COMPOUNDS
        plot_compounds = sorted([c for c in self.compounds if c.lower() != 'glucose'])
        if not plot_compounds:
            return None
            
        bar_data = {c: {'means': [], 'stds': [], 'raw_yields': []} for c in plot_compounds}
        glc_col = self.get_amount_column('Glucose')
        start_glucose = 9.09
        
        for m in mutants:
            key = f"{medium}_{m}"
            if key in medium_mutants:
                raw_df = medium_mutants[key]['raw_data']
                for comp in plot_compounds:
                    comp_col = self.get_amount_column(comp)
                    if comp_col in raw_df.columns and glc_col in raw_df.columns:
                        consumed_glc = start_glucose - raw_df[glc_col]
                        yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                        valid_yields = yields[~np.isnan(yields)]
                        
                        bar_data[comp]['means'].append(np.nanmean(yields) if len(valid_yields) > 0 else 0)
                        bar_data[comp]['stds'].append(np.nanstd(yields) if len(valid_yields) > 0 else 0)
                        bar_data[comp]['raw_yields'].append(valid_yields)
                    else:
                        bar_data[comp]['means'].append(0)
                        bar_data[comp]['stds'].append(0)
                        bar_data[comp]['raw_yields'].append(np.array([]))
            else:
                for comp in plot_compounds:
                    bar_data[comp]['means'].append(0)
                    bar_data[comp]['stds'].append(0)
                    bar_data[comp]['raw_yields'].append(np.array([]))
                    
        fig, ax = plt.subplots(figsize=(14, 7)) 
        x = np.arange(len(mutants))
        width = 0.8 / len(plot_compounds)
        colors = plt.cm.Set2(np.linspace(0, 1, len(plot_compounds)))
        
        max_bar_height = 0
        bars_dict = {}
        
        for i, comp in enumerate(plot_compounds):
            offset = (i - len(plot_compounds)/2 + 0.5) * width
            bars = ax.bar(x + offset, bar_data[comp]['means'], width, yerr=bar_data[comp]['stds'], 
                   label=comp, capsize=3, color=colors[i], edgecolor='black')
            bars_dict[comp] = offset
                 
            for bar, mean_val, std_val in zip(bars, bar_data[comp]['means'], bar_data[comp]['stds']):
                if pd.notna(mean_val) and pd.notna(std_val):
                    # Max Höhe aktualisieren für Signifikanzklammern
                    max_bar_height = max(max_bar_height, mean_val + std_val)
                    
                    # Text über den Balken schreiben (3 Nachkommastellen)
                    ax.annotate(f'{mean_val:.3f}\n±{std_val:.3f}',
                               xy=(bar.get_x() + bar.get_width() / 2, mean_val + std_val),
                               ha='center', va='bottom', 
                               fontsize=self.plot_settings['font_size']-2)
                   
        ax.set_xticks(x)
        
        # Mathematische Formatierung
        display_labels = [r'$\Delta\mathit{' + m + '}$' for m in mutants]
        ax.set_xticklabels(display_labels, rotation=45, ha='right', fontsize=self.plot_settings['font_size'] + 2)
        ax.set_ylabel(f'Yield $Y_{{p/s}}$ [g/g]', fontsize=self.plot_settings['font_size'])
        
        ax.legend(title='Compounds', loc='best', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        
        # Signifikanz-Klammern (Ref: dns)
        if 'dns' in mutants:
            ctrl_idx = mutants.index('dns') 
            y_offset = max_bar_height * 0.15
            current_h = max_bar_height + y_offset
            
            for comp in plot_compounds:
                ctrl_data = bar_data[comp]['raw_yields'][ctrl_idx]
                offset = bars_dict[comp]
                
                if len(ctrl_data) >= 2:
                    for i, m in enumerate(mutants):
                        if i == ctrl_idx: continue
                        
                        test_data = bar_data[comp]['raw_yields'][i]
                        if len(test_data) >= 2:
                            t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                            
                            if p_val < 0.05:
                                text = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
                                x1 = ctrl_idx + offset
                                x2 = i + offset
                                h_tick = max_bar_height * 0.02
                                
                                ax.plot([x1, x1, x2, x2], [current_h, current_h+h_tick, current_h+h_tick, current_h], lw=1.2, c='black')
                                ax.text((x1+x2)*.5, current_h+h_tick, text, ha='center', va='bottom', color='black', fontsize=self.plot_settings['font_size'] + 2)
                                current_h += (max_bar_height * 0.08)
                                
            ax.set_ylim(0, current_h + (max_bar_height * 0.05))
        
        plt.tight_layout()
        filepath = str(Path(output_dir) / f"biolector_{medium}_all_compounds_yield.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'], bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        print(f"✓ Yield-Plot (alle Compounds, signifikant) erstellt: {filepath}")
        return filepath

    def plot_biolector_yield_all_compounds_media(self, biolector_data: Dict, mutant: str, output_dir: str = './plots') -> Optional[str]:
        """Plottet alle Compounds (Yield) über alle Medien für eine bestimmte Mutante inkl. Signifikanz (Ref: M1)."""
        if not biolector_data:
            return None
            
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        media_order = ['M1', 'FeCN', 'ABTS']
        display_labels = ['M1', '5 mM FeCN', '100 \u03BCM ABTS']
        
        # Prüfen, ob die Mutante überhaupt existiert
        mutant_exists = any(v.get('mutant') == mutant for v in biolector_data.values())
        if not mutant_exists:
            return None
            
        # ALPHABETISCHE SORTIERUNG COMPOUNDS
        plot_compounds = sorted([c for c in self.compounds if c.lower() != 'glucose'])
        if not plot_compounds:
            return None
            
        bar_data = {c: {'means': [], 'stds': [], 'raw_yields': []} for c in plot_compounds}
        glc_col = self.get_amount_column('Glucose')
        start_glucose = 9.09
        
        for medium in media_order:
            key = f"{medium}_{mutant}"
            if key in biolector_data:
                raw_df = biolector_data[key]['raw_data']
                for comp in plot_compounds:
                    comp_col = self.get_amount_column(comp)
                    if comp_col in raw_df.columns and glc_col in raw_df.columns:
                        consumed_glc = start_glucose - raw_df[glc_col]
                        yields = np.where(consumed_glc > 0, raw_df[comp_col] / consumed_glc, np.nan)
                        valid_yields = yields[~np.isnan(yields)]
                        
                        bar_data[comp]['means'].append(np.nanmean(yields) if len(valid_yields) > 0 else 0)
                        bar_data[comp]['stds'].append(np.nanstd(yields) if len(valid_yields) > 0 else 0)
                        bar_data[comp]['raw_yields'].append(valid_yields)
                    else:
                        bar_data[comp]['means'].append(0)
                        bar_data[comp]['stds'].append(0)
                        bar_data[comp]['raw_yields'].append(np.array([]))
            else:
                for comp in plot_compounds:
                    bar_data[comp]['means'].append(0)
                    bar_data[comp]['stds'].append(0)
                    bar_data[comp]['raw_yields'].append(np.array([]))
                    
        fig, ax = plt.subplots(figsize=(14, 7)) 
        x = np.arange(len(media_order))
        width = 0.8 / len(plot_compounds)
        colors = plt.cm.Set2(np.linspace(0, 1, len(plot_compounds)))
        
        max_bar_height = 0
        bars_dict = {}
        
        for i, comp in enumerate(plot_compounds):
            offset = (i - len(plot_compounds)/2 + 0.5) * width
            bars = ax.bar(x + offset, bar_data[comp]['means'], width, yerr=bar_data[comp]['stds'], 
                   label=comp, capsize=3, color=colors[i], edgecolor='black')
            bars_dict[comp] = offset
            
            for bar, mean_val, std_val in zip(bars, bar_data[comp]['means'], bar_data[comp]['stds']):
                if pd.notna(mean_val) and pd.notna(std_val):
                    max_bar_height = max(max_bar_height, mean_val + std_val)
                    
                    ax.annotate(f'{mean_val:.3f}\n±{std_val:.3f}',
                               xy=(bar.get_x() + bar.get_width() / 2, mean_val + std_val),
                               ha='center', va='bottom', 
                               fontsize=self.plot_settings['font_size']-2)
                   
        ax.set_xticks(x)
        ax.set_xticklabels(display_labels, rotation=0, fontsize=self.plot_settings['font_size'] + 2)
        ax.set_ylabel(f'Yield $Y_{{p/s}}$ [g/g]', fontsize=self.plot_settings['font_size'])
        
        # Legende INNERHALB des Diagramms (loc='best' sucht den freiesten Platz)
        ax.legend(title='Compounds', loc='best', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        
        # =========================================================
        # SIGNIFIKANZ-KLAMMERN EINZEICHNEN (Ref: M1)
        # =========================================================
        ctrl_idx = 0 # M1 ist immer Index 0
        y_offset = max_bar_height * 0.15
        current_h = max_bar_height + y_offset
        
        for comp in plot_compounds:
            ctrl_data = bar_data[comp]['raw_yields'][ctrl_idx]
            offset = bars_dict[comp]
            
            if len(ctrl_data) >= 2:
                for i, medium in enumerate(media_order):
                    if i == ctrl_idx: continue
                    
                    test_data = bar_data[comp]['raw_yields'][i]
                    if len(test_data) >= 2:
                        t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                        
                        if p_val < 0.05:
                            text = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
                            x1 = ctrl_idx + offset
                            x2 = i + offset
                            h_tick = max_bar_height * 0.02
                            
                            ax.plot([x1, x1, x2, x2], [current_h, current_h+h_tick, current_h+h_tick, current_h], lw=1.2, c='black')
                            ax.text((x1+x2)*.5, current_h+h_tick, text, ha='center', va='bottom', color='black', fontsize=self.plot_settings['font_size'] + 2)
                            
                            current_h += (max_bar_height * 0.08)
                            
        ax.set_ylim(0, current_h + (max_bar_height * 0.05))
        
        plt.tight_layout()
        filepath = str(Path(output_dir) / f"biolector_yield_all_compounds_media_{mutant}.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'], bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        print(f"✓ Yield-Plot (alle Compounds, alle Medien, signifikant) erstellt: {filepath}")
        return filepath

    def plot_biolector_medium_comparison_yield(self, biolector_data: Dict, compound: str, mutant: str, output_dir: str = './plots') -> Optional[str]:
        """Vergleicht die Medien für eine Mutante bezüglich des Yields (g/g) inkl. Signifikanz (Ref: M1)."""
        if not biolector_data:
            return None
        if compound.lower() == 'glucose':
            return None  
            
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        media_order = ['M1', 'FeCN', 'ABTS']
        display_labels = ['M1', '5 mM FeCN', '100 \u03BCM ABTS']
        colors = ['#0072B2', '#E69F00', '#CC79A7']
        
        yield_means = []
        yield_stds = []
        raw_yields_dict = {}
        
        amount_col = self.get_amount_column(compound)
        glc_col = self.get_amount_column('Glucose')
        start_glucose = 9.09
        
        for medium in media_order:
            key = f"{medium}_{mutant}"
            if key in biolector_data:
                raw_df = biolector_data[key]['raw_data']
                if amount_col in raw_df.columns and glc_col in raw_df.columns:
                    consumed_glc = start_glucose - raw_df[glc_col]
                    yields = np.where(consumed_glc > 0, raw_df[amount_col] / consumed_glc, np.nan)
                    valid_yields = yields[~np.isnan(yields)]
                    
                    yield_means.append(np.nanmean(valid_yields) if len(valid_yields)>0 else 0)
                    yield_stds.append(np.nanstd(valid_yields) if len(valid_yields)>0 else 0)
                    raw_yields_dict[medium] = valid_yields
                else:
                    yield_means.append(0)
                    yield_stds.append(0)
                    raw_yields_dict[medium] = np.array([])
            else:
                yield_means.append(0)
                yield_stds.append(0)
                raw_yields_dict[medium] = np.array([])
                
        if all(v == 0 for v in yield_means):
            return None
            
        fig, ax = plt.subplots(figsize=self.plot_settings['figsize'])
        x_pos = np.arange(len(media_order))
        
        bars = ax.bar(x_pos, yield_means, yerr=yield_stds, capsize=5,
                     color=colors, edgecolor='black', linewidth=1.5, alpha=0.8)
                     
        for bar, val, std in zip(bars, yield_means, yield_stds):
            if pd.notna(val) and pd.notna(std):
                # Text platzieren (y-Position = Balkenhöhe + Standardabweichung)
                ax.annotate(f'{val:.3f}\n±{std:.3f}',
                           xy=(bar.get_x() + bar.get_width() / 2, val + std),
                           ha='center', va='bottom', 
                           fontsize=self.plot_settings['font_size']-2)
                           
        ax.set_xticks(x_pos)
        ax.set_xticklabels(display_labels, fontsize=self.plot_settings['font_size'])
        ax.set_ylabel(f'{compound} Yield $Y_{{p/s}}$ [g/g]', fontsize=self.plot_settings['font_size'])
        ax.grid(True, alpha=0.3, axis='y')
        
        max_bar_height = max([m + s for m, s in zip(yield_means, yield_stds) if pd.notna(m) and pd.notna(s)] + [0])
        
        # Signifikanztest gegen 'M1'
        ctrl_idx = 0
        ctrl_data = raw_yields_dict['M1']
        
        if len(ctrl_data) >= 2:
            y_offset = max_bar_height * 0.15
            current_h = max_bar_height + y_offset
            
            for i, medium in enumerate(media_order):
                if i == ctrl_idx: continue
                test_data = raw_yields_dict[medium]
                if len(test_data) >= 2:
                    t_stat, p_val = stats.ttest_ind(ctrl_data, test_data, equal_var=False)
                    if p_val < 0.05:
                        text = '***' if p_val < 0.001 else '**' if p_val < 0.01 else '*'
                        x1, x2 = ctrl_idx, i
                        h_tick = max_bar_height * 0.02
                        ax.plot([x1, x1, x2, x2], [current_h, current_h+h_tick, current_h+h_tick, current_h], lw=1.2, c='black')
                        ax.text((x1+x2)*.5, current_h+h_tick, text, ha='center', va='bottom', color='black', fontsize=self.plot_settings['font_size'] + 2)
                        current_h += (max_bar_height * 0.12)
            ax.set_ylim(0, current_h + (max_bar_height * 0.05))

        plt.tight_layout()
        filepath = str(Path(output_dir) / f"biolector_yield_{mutant}_{compound}_media_comparison.png")
        plt.savefig(filepath, dpi=self.plot_settings['dpi'], bbox_inches='tight', format=self.plot_settings['save_format'])
        plt.close()
        print(f"✓ Yield-Medien-Vergleich Plot erstellt: {filepath}")
        return filepath
    
    def export_all_tables(self, output_dir: str = './tables') -> Dict[str, str]:
        """Exportiert alle Tabellen als CSV"""
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        exported = {}
        
        ecuev = self.group_ecuev_data()
        if ecuev['mutants_mean']:
            df = self.create_summary_table(ecuev['mutants_mean'], 'Ecuev Mutanten')
            filepath = f"{output_dir}/ecuev_mutants_summary.csv"
            df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
            exported['ecuev_mutants'] = filepath
            print(f"✓ Tabelle: {filepath}")
            
            ecuev_yield_df = self.calculate_ecuev_yields(ecuev)
            if ecuev_yield_df is not None:
                yield_filepath = f"{output_dir}/ecuev_yields_summary.csv"
                ecuev_yield_df.to_csv(yield_filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
                exported['ecuev_yields'] = yield_filepath
                print(f"✓ Yield-Tabelle: {yield_filepath}")
                
            ecuev_stats_path = self.export_ecuev_statistics(ecuev, output_dir)
            if ecuev_stats_path:
                exported['ecuev_stats'] = ecuev_stats_path
                print(f"✓ Statistik-Tabelle (Ecuev): {ecuev_stats_path}")
        
        if ecuev['controls']:
            df = self.create_summary_table(ecuev['controls'], 'Ecuev Kontrollen')
            filepath = f"{output_dir}/ecuev_controls_summary.csv"
            df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
            exported['ecuev_controls'] = filepath
            print(f"✓ Tabelle: {filepath}")
        
        dated = self.group_dated_data()
        if dated:
            # 1. Normale Konzentrationstabelle exportieren
            df = self.create_summary_table(dated, 'Datierte Experimente')
            filepath = f"{output_dir}/dated_experiments_summary.csv"
            df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
            exported['dated'] = filepath
            print(f"✓ Tabelle: {filepath}")
            
            # 2. Spezifische Yield-Tabelle exportieren
            yield_df = self.calculate_dated_yields()
            if yield_df is not None:
                yield_filepath = f"{output_dir}/dated_yields_summary.csv"
                yield_df_rounded = yield_df.copy()
                numeric_cols = yield_df_rounded.select_dtypes(include=[np.number]).columns
                yield_df_rounded[numeric_cols] = yield_df_rounded[numeric_cols].round(4)
                yield_df_rounded.to_csv(yield_filepath, index=False, sep=';', decimal=',')
                exported['dated_yields'] = yield_filepath
                print(f"✓ Yield-Tabelle: {yield_filepath}")
                
                
                # Zusaetzlich: Yield-Zeitverlauf mit jedem einzelnen Messzeitpunkt exportieren
                yield_time_df = self.calculate_dated_yield_timecourse()
                if yield_time_df is not None:
                    yield_time_path = f"{output_dir}/dated_yields_timecourse.csv"
                    yield_time_export = yield_time_df.copy()
                    numeric_cols = yield_time_export.select_dtypes(include=[np.number]).columns
                    yield_time_export[numeric_cols] = yield_time_export[numeric_cols].round(4)
                    yield_time_export.to_csv(yield_time_path, index=False, sep=';', decimal=',', encoding='utf-8-sig')
                    exported['dated_yields_timecourse'] = yield_time_path
                    print(f"✓ Yield-Zeitverlauf-Tabelle: {yield_time_path}")

                stats_filepath = self.export_yield_statistics(yield_df, output_dir)
                if stats_filepath:
                    exported['dated_yields_stats'] = stats_filepath
                    print(f"✓ Statistik-Tabelle (p-Werte): {stats_filepath}")

            # 3. Endkonzentrationen am letzten Zeitpunkt vergleichen (M1 vs FeCN/ABTS)
            end_conc_stats_path = self.export_dated_end_concentration_statistics(output_dir)
            if end_conc_stats_path:
                exported['dated_end_concentration_stats'] = end_conc_stats_path
                print(f"✓ Endkonzentrations-Statistik: {end_conc_stats_path}")
        
        biolector = self.group_biolector_data()
        if biolector:
            df = self.create_summary_table(biolector, 'Biolector Experimente')
            filepath = f"{output_dir}/biolector_summary.csv"
            df.to_csv(filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
            exported['biolector'] = filepath
            print(f"✓ Tabelle: {filepath}")
            
            biolector_yield_df = self.calculate_biolector_yields(biolector)
            if biolector_yield_df is not None:
                yield_filepath = f"{output_dir}/biolector_yields_summary.csv"
                biolector_yield_df.to_csv(yield_filepath, index=False, sep=';', decimal=',', encoding='utf-8-sig')
                exported['biolector_yields'] = yield_filepath
                print(f"✓ Yield-Tabelle: {yield_filepath}")
        
            stats_filepath = self.export_biolector_statistics(biolector, output_dir)
            if stats_filepath:
                exported['biolector_stats'] = stats_filepath
                print(f"✓ Statistik-Tabelle (p-Werte): {stats_filepath}")
        
        print(f"\n✓ {len(exported)} Tabellen exportiert")
        return exported


class HPLCAnalyzerGUI:
    """Graphische Benutzeroberfläche"""
    
    def __init__(self, root):
        self.root = root
        self.root.title("HPLC Data Analyzer - Universelles Tool")
        self.root.geometry("1200x800")
        
        self.analyzer = HPLCDataAnalyzer()
        
        self.setup_ui()
    
    def setup_ui(self):
        """UI aufbauen"""
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        # File Selection
        file_frame = ttk.LabelFrame(main_frame, text="1. Datei auswählen", padding="10")
        file_frame.grid(row=0, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.file_path_var = tk.StringVar()
        ttk.Entry(file_frame, textvariable=self.file_path_var, width=80).grid(row=0, column=0, padx=5)
        ttk.Button(file_frame, text="Durchsuchen...", command=self.browse_file).grid(row=0, column=1)
        ttk.Button(file_frame, text="Laden", command=self.load_file).grid(row=0, column=2, padx=10)
        
        # Compound Selection
        compound_frame = ttk.LabelFrame(main_frame, text="2. Stoffe konfigurieren", padding="10")
        compound_frame.grid(row=1, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.compound_listbox = tk.Listbox(compound_frame, height=6, selectmode=tk.MULTIPLE)
        self.compound_listbox.grid(row=0, column=0, rowspan=4, sticky=(tk.W, tk.E))
        
        for compound in self.analyzer.compounds:
            self.compound_listbox.insert(tk.END, compound)
        
        btn_frame = ttk.Frame(compound_frame)
        btn_frame.grid(row=0, column=1, sticky=tk.N)
        
        ttk.Button(btn_frame, text="Hinzufügen", command=self.add_compound).pack(pady=2)
        ttk.Button(btn_frame, text="Entfernen", command=self.remove_compound).pack(pady=2)
        
        self.new_compound_var = tk.StringVar()
        ttk.Entry(btn_frame, textvariable=self.new_compound_var, width=15).pack(pady=2)
        
        # Analysis Options
        analysis_frame = ttk.LabelFrame(main_frame, text="3. Analyse-Optionen", padding="10")
        analysis_frame.grid(row=2, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.var_ecuev = tk.BooleanVar(value=True)
        self.var_dated = tk.BooleanVar(value=True)
        self.var_biolector = tk.BooleanVar(value=True)
        
        ttk.Checkbutton(analysis_frame, text="Ecuev-Experimente", variable=self.var_ecuev).grid(row=0, column=0, padx=10)
        ttk.Checkbutton(analysis_frame, text="Datierte Experimente (260707-260714)", variable=self.var_dated).grid(row=0, column=1, padx=10)
        ttk.Checkbutton(analysis_frame, text="Biolector-Experimente", variable=self.var_biolector).grid(row=0, column=2, padx=10)
        
        # Plot Options
        plot_frame = ttk.LabelFrame(main_frame, text="4. Diagramm-Optionen", padding="10")
        plot_frame.grid(row=3, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.var_timecourse = tk.BooleanVar(value=True)
        self.var_medium_comp = tk.BooleanVar(value=True)
        self.var_mutant_comp = tk.BooleanVar(value=True)
        
        ttk.Checkbutton(plot_frame, text="Zeitverlauf pro Stoff", variable=self.var_timecourse).grid(row=0, column=0, padx=10)
        ttk.Checkbutton(plot_frame, text="Medien-Vergleich", variable=self.var_medium_comp).grid(row=0, column=1, padx=10)
        ttk.Checkbutton(plot_frame, text="Mutanten-Vergleich (Biolector)", variable=self.var_mutant_comp).grid(row=0, column=2, padx=10)
        
        # Output Directory
        output_frame = ttk.LabelFrame(main_frame, text="5. Ausgabe-Verzeichnis", padding="10")
        output_frame.grid(row=4, column=0, columnspan=3, sticky=(tk.W, tk.E), pady=5)
        
        self.output_dir_var = tk.StringVar(value='./output')
        ttk.Entry(output_frame, textvariable=self.output_dir_var, width=80).grid(row=0, column=0, padx=5)
        ttk.Button(output_frame, text="Ändern", command=self.browse_output).grid(row=0, column=1)
        
        # Run Button
        run_frame = ttk.Frame(main_frame)
        run_frame.grid(row=5, column=0, columnspan=3, pady=20)
        
        ttk.Button(run_frame, text="ANALYSE STARTEN", command=self.run_analysis,
                  style='Accent.TButton').pack(side=tk.LEFT, padx=10)
        ttk.Button(run_frame, text="Vorschau", command=self.preview_data).pack(side=tk.LEFT, padx=10)
        
        # Status/Log
        log_frame = ttk.LabelFrame(main_frame, text="Status-Log", padding="10")
        log_frame.grid(row=6, column=0, columnspan=3, sticky=(tk.W, tk.E, tk.N, tk.S), pady=5)
        
        self.log_text = tk.Text(log_frame, height=15, width=120)
        self.log_text.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
        
        scrollbar = ttk.Scrollbar(log_frame, orient=tk.VERTICAL, command=self.log_text.yview)
        scrollbar.grid(row=0, column=1, sticky=(tk.N, tk.S))
        self.log_text.configure(yscrollcommand=scrollbar.set)
        
        self.root.columnconfigure(0, weight=1)
        self.root.rowconfigure(0, weight=1)
        main_frame.columnconfigure(0, weight=1)
        main_frame.rowconfigure(6, weight=1)
    
    def log(self, message: str):
        timestamp = datetime.now().strftime("%H:%M:%S")
        self.log_text.insert(tk.END, f"[{timestamp}] {message}\n")
        self.log_text.see(tk.END)
    
    def browse_file(self):
        filename = filedialog.askopenfilename(
            title="HPLC Excel-Datei auswählen",
            # Hinzugefügt: CSV-Unterstützung im Auswahlmenü
            filetypes=[
                ("Excel-Dateien", "*.xlsx *.xls"),
                ("Excel 2007-365", "*.xlsx"),
                ("Excel 97-2003", "*.xls"),
                ("CSV-Dateien", "*.csv"),
                ("Alle Dateien", "*.*")
            ]
        )
        if filename:
            self.file_path_var.set(filename)
    
    def browse_output(self):
        directory = filedialog.askdirectory(title="Ausgabeverzeichnis wählen")
        if directory:
            self.output_dir_var.set(directory)
    
    def add_compound(self):
        compound = self.new_compound_var.get().strip()
        if compound:
            if self.analyzer.add_compound(compound):
                self.compound_listbox.insert(tk.END, compound)
                self.new_compound_var.set("")
                self.log(f"✓ Stoff hinzugefügt: {compound}")
    
    def remove_compound(self):
        selection = self.compound_listbox.curselection()
        if selection:
            for index in reversed(selection):
                compound = self.compound_listbox.get(index)
                self.analyzer.remove_compound(compound)
                self.compound_listbox.delete(index)
                self.log(f"✓ Stoff entfernt: {compound}")
    
    def load_file(self):
        filepath = self.file_path_var.get()
        if not filepath:
            messagebox.showerror("Fehler", "Bitte Datei auswählen!")
            return
        
        if self.analyzer.load_data(filepath):
            self.log(f"✓ Datei geladen: {filepath}")
            self.log(f"✓ {len(self.analyzer.data)} Zeilen, {len(self.analyzer.data.columns)} Spalten")
        else:
            # Klarere Fehlermeldung
            messagebox.showerror("Fehler", "Datei konnte nicht geladen werden!\nStelle sicher, dass die Datei nicht in Excel geöffnet ist.")
    
    def preview_data(self):
        if self.analyzer.data is None:
            messagebox.showinfo("Info", "Bitte zuerst Datei laden!")
            return
        
        preview_win = tk.Toplevel(self.root)
        preview_win.title("Datenvorschau")
        preview_win.geometry("800x600")
        
        text = tk.Text(preview_win)
        text.pack(fill=tk.BOTH, expand=True)
        
        text.insert(tk.END, self.analyzer.data.head(20).to_string())
    
    def run_analysis(self):
        if self.analyzer.data is None:
            messagebox.showerror("Fehler", "Bitte zuerst Datei laden!")
            return
        
        self.log("=" * 60)
        self.log("ANALYSE GESTARTET")
        self.log("=" * 60)
        
        output_dir = self.output_dir_var.get()
        
        try:
            self.log("\n📊 Tabellen werden erstellt...")
            tables = self.analyzer.export_all_tables(str(Path(output_dir) / "tables"))
            for name, path in tables.items():
                self.log(f"  ✓ {name}: {path}")
            
            self.log("\n📈 Diagramme werden erstellt...")
            
            dated_data = self.analyzer.group_dated_data() if self.var_dated.get() else {}
            biolector_data = self.analyzer.group_biolector_data() if self.var_biolector.get() else {}
            
            plots_dir = Path(output_dir) / "plots"
            plots_dir.mkdir(parents=True, exist_ok=True)
            
            # Timed Data: Alle Compounds gemeinsam pro Medium
            # Nur EINMAL pro Medium erzeugen, nicht innerhalb der Compound-Schleife.
            if self.var_timecourse.get() and dated_data:
                self.log("\n  Erstelle Gesamt-Zeitverläufe aller Compounds...")
                for medium in ['M1', 'FeCN', 'ABTS']:
                    self.analyzer.plot_all_compounds_by_medium(
                        dated_data, medium, str(plots_dir))

            # Alphabetische Sortierung der Compounds für die Diagrammerstellung erzwungen
            for compound in sorted(self.analyzer.compounds):
                self.log(f"\n  Stoff: {compound}")

                if self.var_timecourse.get() and dated_data:
                    self.analyzer.plot_compound_timecourse(
                        dated_data, compound, str(plots_dir))

                    # Für alle Produkte außer Glucose zusätzlich den kumulativen Yield-Zeitverlauf plotten.
                    if compound.lower() != 'glucose':
                        yield_time_df = self.analyzer.calculate_dated_yield_timecourse()
                        if yield_time_df is not None:
                            self.analyzer.plot_dated_yield_timecourse(
                                yield_time_df, compound, str(plots_dir))

                if self.var_medium_comp.get() and dated_data:
                    self.analyzer.plot_medium_comparison(
                        dated_data, compound, str(plots_dir))

                if self.var_mutant_comp.get() and biolector_data:
                    # 1. Plot: Konzentration [g/L] für Mutanten in einem Medium vergleichen
                    for medium in ['M1', 'FeCN', 'ABTS']:
                        self.analyzer.plot_biolector_mutant_comparison(
                            biolector_data, compound, medium, str(plots_dir))

                    # 2. Plot: Yield [g/g] über verschiedene Medien für alle Mutanten vergleichen
                    mutants = list(set([
                        v['mutant'] for v in biolector_data.values() if v.get('mutant')
                    ]))
                    for mutant in mutants:
                        self.analyzer.plot_biolector_medium_comparison_yield(
                            biolector_data, compound, mutant, str(plots_dir))

            # 3. Plot: Yields (g/g) ALLER Compounds zusammengefasst pro Medium
            if self.var_mutant_comp.get() and biolector_data:
                self.log("\n  Erstelle BioLector Yield-Plots (g/g) pro Medium...")
                for medium in ['M1', 'FeCN', 'ABTS']:
                    self.analyzer.plot_biolector_yield_all_compounds(
                        biolector_data, medium, str(plots_dir))

            # 4. Plot: Yields (g/g) ALLER Compounds zusammengefasst pro Mutante (Medien-Vergleich)
            if self.var_mutant_comp.get() and biolector_data:
                self.log("\n  Erstelle BioLector Yield-Plots (g/g) als Medien-Vergleich...")
                mutants = list(set([
                    v['mutant'] for v in biolector_data.values() if v.get('mutant')
                ]))
                for mutant in mutants:
                    self.analyzer.plot_biolector_yield_all_compounds_media(
                        biolector_data, mutant, str(plots_dir))

            # Yield-Plots generieren
            if self.var_medium_comp.get() and dated_data:
                self.log("\n  Erstelle berechnete Yield-Plots (g/g)...")
                yield_df = self.analyzer.calculate_dated_yields()
                if yield_df is not None:
                    yield_plots = self.analyzer.plot_dated_yields(
                        yield_df, str(plots_dir))
                    for p in yield_plots:
                        self.log(f"  ✓ {Path(p).name}")

                    combined_plot = self.analyzer.plot_dated_yields_m1_fecn_all_compounds(
                        yield_df, str(plots_dir))
                    if combined_plot:
                        self.log(f"  ✓ {Path(combined_plot).name}")

            self.log("\n" + "=" * 60)
            self.log("✅ ANALYSE ERFOLGREICH ABGESCHLOSSEN")
            self.log("=" * 60)
            messagebox.showinfo("Erfolg", f"Analyse abgeschlossen!\nAusgabe: {output_dir}")
            
        except Exception as e:
            self.log(f"\n✗ FEHLER: {e}")
            import traceback
            traceback.print_exc()
            messagebox.showerror("Fehler", f"Analyse fehlgeschlagen:\n{e}")

def main():
    """Hauptfunktion"""
    print("=" * 60)
    print("HPLC Data Analyzer - Start")
    print("=" * 60)
    
    root = tk.Tk()
    app = HPLCAnalyzerGUI(root)
    root.mainloop()

if __name__ == "__main__":
    main()
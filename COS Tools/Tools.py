from glob import glob
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from datetime import datetime, timedelta
from gspread_formatting import set_frozen
from openpyxl import load_workbook
from numpy import arcsin, sin, cos, sqrt
from math import pi
import os
import time
import gspread
import pandas as pd
import requests as rq
import urllib3
import queue
import threading
import sys
import ctypes


# Disable SSL warnings
urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

# Define ctypes structures and constants
LPSECURITY_ATTRIBUTES = ctypes.wintypes.LPVOID
HANDLE = ctypes.wintypes.HANDLE
LPCTSTR = ctypes.wintypes.LPCWSTR
LPTSTR = ctypes.wintypes.LPWSTR
DWORD = ctypes.wintypes.DWORD

class SECURITY_ATTRIBUTES(ctypes.Structure):
    _fields_ = [("nLength", DWORD),
                ("lpSecurityDescriptor", LPSECURITY_ATTRIBUTES),
                ("bInheritHandle", ctypes.wintypes.BOOL)]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.nLength = ctypes.sizeof(self)

# Windows API functions
CreateMutex = ctypes.windll.kernel32.CreateMutexW
CreateMutex.argtypes = [LPSECURITY_ATTRIBUTES, ctypes.wintypes.BOOL, LPCTSTR]
CreateMutex.restype = HANDLE

CloseHandle = ctypes.windll.kernel32.CloseHandle
CloseHandle.argtypes = [HANDLE]
CloseHandle.restype = ctypes.wintypes.BOOL

GetLastError = ctypes.windll.kernel32.GetLastError
GetLastError.argtypes = []
GetLastError.restype = DWORD

MessageBox = ctypes.windll.user32.MessageBoxW
MessageBox.argtypes = [ctypes.wintypes.HWND, LPCTSTR, LPCTSTR, ctypes.wintypes.UINT]
MessageBox.restype = ctypes.wintypes.INT


# Tooltip
class Tooltip:
    def __init__(self, widget, text):
        self.widget = widget
        self.text = text
        self.tooltip_window = None

    def show_tooltip(self, event):
        x, y, _, _ = self.widget.bbox(tk.CURRENT)
        x_root, y_root = self.widget.winfo_rootx(), self.widget.winfo_rooty()
        self.tooltip_window = tk.Toplevel(self.widget)
        self.tooltip_window.wm_overrideredirect(True)
        self.tooltip_window.wm_geometry(f"+{x_root + x + 10}+{y_root + y + 10}")
        label = tk.Label(self.tooltip_window, text=self.text, background="lightyellow", relief="solid", borderwidth=1)
        label.pack()

    def hide_tooltip(self, event):
        if self.tooltip_window:
            self.tooltip_window.destroy()
            self.tooltip_window = None


# Calculates the distance of a station to an epicenter
def st_dist(df, x, y):
    dist = 6371*2*arcsin(sqrt(sin((df['Latitude']*pi/180 - x*pi/180)/2)**2 + 
                        cos(x*pi/180)*cos(df['Latitude']*pi/180)*sin((df['Longitude']*pi/180 - y*pi/180)/2)**2))
    
    return dist


# Function to compare the list of stations used to list of related stations to the eq
def compare_stations_used_to_quake(stations_used, stations_related_to_quake, filename):

    # Find the first occurrence of each station in stations_related_to_quake
    quake_indices = []
    incorrect_stations = []
    stations_with_asterisk = []
    stations_vol = []
    stations_not_found = []
    for station in stations_used:
        try:
            quake_indices.append(stations_related_to_quake.index(station[:3]))
            
            if station.endswith('*'):
                stations_with_asterisk.append(station[:3])
                
            if station.startswith('V'):
                stations_vol.append(station[:3])
                

        except ValueError:
            station = station.rstrip('*')
            stations_not_found.append((f"{station}", f"{filename}"))
            
    # Index for vol stations
    vols = 0
    # Check if the indices are in ascending order in the quake_indices list
    for i in range(len(quake_indices) - 1):
        
        if i == 0:
            if quake_indices[i] >= quake_indices[i + 1] and abs(quake_indices[i] - quake_indices[i + 1]) != 1 and (stations_related_to_quake[quake_indices[i]][:2] if stations_related_to_quake[quake_indices[i]].startswith('V') else stations_related_to_quake[quake_indices[i]]) in stations_with_asterisk:
                incorrect_stations.append(stations_related_to_quake[quake_indices[i]])
        else:
    
            if quake_indices[i] > quake_indices[i + 1] and abs(quake_indices[i] - quake_indices[i + 1]) != 1:

                if quake_indices[i + 1] > quake_indices[i-1] and abs(quake_indices[i + 1] - quake_indices[i-1] != 1) and stations_related_to_quake[quake_indices[i]] in stations_with_asterisk:
                    incorrect_stations.append(stations_vol[vols] if stations_related_to_quake[quake_indices[i]].startswith('V') else stations_related_to_quake[quake_indices[i]])

                
                elif (quake_indices[i-1] < quake_indices[i] 
                        and stations_related_to_quake[quake_indices[i+1]] in stations_with_asterisk 
                        or stations_related_to_quake[quake_indices[i]] in incorrect_stations
                    ):
                    
                    incorrect_stations.append(stations_vol[vols] if stations_related_to_quake[quake_indices[i+1]].startswith('V') else stations_related_to_quake[quake_indices[i+1]])
              
        if stations_related_to_quake[quake_indices[i]].startswith('V'):
            vols += 1
    
   
    return incorrect_stations, stations_not_found


# Function to check the file
def check_eqp_file(file_path):
        with open(file_path, 'r') as eqp_file:
            pri_values = []       # List to store values from the "pri" column
            sec_values = []       # List to store values from the "sec" column
            sta_unique = set()
            sta_duplicate = []    # Set to store unique values from the "sta" column
            sta_check = True      # Check for duplicates
            sta_used = []
            
            sta_vol = [] # List for volcano
            
            resp_values = [] # List for RESp > 1
            ress_values = [] # List for RESs > 1
            failed_sta_resp = [] # List for failed station RESp
            failed_sta_ress = [] # List for failed station RESs

            # Variables to store values
            rms_value = None
            mag_value = None
            seconds_value = None
            depth_value = None
            resp_value = None
            ress_value = None
            
            # Initialize sections
            data_section = False  # Flag to indicate if we're inside the [Data] section
            rms_section = False
            result_section = False
            res_section = False
            mag_section = False
            issuedby_section = False

            for line in eqp_file:
                line = line.strip()  # Remove leading/trailing whitespaces

                # Check if we're inside the [Data] section
                if line.startswith("[Data]"):
                    data_section = True
                    continue
                
                # Inside the [Result] section, look for lines with RMS, SEC, RESp, and RESs values
                if line.startswith("[Result]"):
                    data_section = False
                    continue
                
                if line.startswith("IT"):
                    rms_section = True
                    continue
                
                # Inside the [Result] section, look for lines with RMS, SEC, RESp, and RESs values
                if line.startswith("ERR"):
                    result_section = True
                    rms_section = False
                    continue
                
                # Inside the STA section, look for lines with RMS, SEC, RESp, and RESs values
                if line.startswith("STA"):
                    res_section = True
                    continue
                
                # Inside the [Magnitude] section, look for lines with RMS, SEC, RESp, and RESs values
                if line.startswith("[Magnitude]"):
                    mag_section = True
                    res_section = False
                    continue
                
                if line.startswith("[IssuedBy]"):
                    issuedby_section = True
                    mag_section = False
                    continue
                
                # Inside the [Data] section, look for lines with pri, sec, and sta values
                if data_section and line:
                    columns = line.split()
                    
                    if len(columns) >= 2:
                        if len(columns) == 2:
                            columns.insert(2, str(0))
                            columns.insert(3, str(0))
                                 
                        elif len(columns) == 3:
                            columns.insert(3, str(0))
                                           
                        sta_value = columns[0]
                        
                        # Skip lines that contain 'sta'
                        if sta_value == 'sta':
                            continue

                        pri_value = columns[1]
                        sec_value = columns[2]
                        
                        try:
                           
                            # Check if pri and sec values meet the criteria
                            if '*' not in pri_value:
                                pri_values.append(float(pri_value))
                            else:
                                sta_used.append(sta_value+'*')
   
                            if '*' not in sec_value and sec_value != '0' and '*' not in pri_value:
                                sec_values.append(float(sec_value))
                                
                            if not any(sta_value == s[:-1] for s in sta_used):
                                sta_used.append(sta_value)
                                
                            
                        except ValueError:
                            pass

                        # Check if the first 3 letters of sta_value are unique within the set
                        
                        sta_prefix = sta_value[:3]
                        if sta_prefix.startswith("V"):
                            sta_vol.append(sta_prefix)
                            
                        else:
                            if sta_prefix in sta_unique:
                                sta_check = False
                                sta_duplicate.append(sta_prefix)
                            else:
                                sta_unique.add(sta_prefix)
                       
    
                if rms_section:
                    values = line.split()
                    
                    if len(values) == 8 :
                        try:
                            rms_value = float(values[7])  # RMS is at index 7
                            depth_value = float(values[3])
                            lat_value = float(values[1])
                            long_value = float(values[2])
                        except ValueError:
                            pass

                    

                if result_section:
                    values = line.split()
                    
                    if len(values) == 4 :
                        try:
                            seconds_value = float(values[3])  # SEC is at index 3
                            result_section = False
                            
                        except ValueError:
                            pass
                        
                   
                if res_section:
                    values = line.split()
                    if len(values) == 10 :
                        try:
                
                            resp_value = float(values[4])  # RESp is at index 4
                            if abs(resp_value) >= 1:
                                resp_values.append(resp_value)
                                failed_sta_resp.append(values[0]) # Append failed station

                            ress_value = float(values[5])  # RESs is at index 5
                            if abs(ress_value) >= 1:
                                ress_values.append(ress_value)
                                failed_sta_ress.append(values[0]) # Append failed station
                            
                        except ValueError:
                            pass

                if mag_section:
                    values = line.split('/')
                    
                    try:
                        
                        mag_value = float(values[0][3:])
                        
                    except ValueError:
                        pass
                    
                    
                if issuedby_section:
                    try:
                        issuedby_value = line
                        
                    except ValueError:
                        pass
                        
            
            sta_vol_counts = {}
            for station in sta_vol:
                prefix = station[:2]
                sta_vol_counts[prefix] = sta_vol_counts.get(prefix, 0) + 1
            
            # Filter vol_stations with counts greater than 2
            sta_vol_duplicates = [station[:2] for station in sta_vol if sta_vol_counts.get(station[:2], 0) > 2]
            
            if sta_vol_duplicates:
                sta_check = False
 
            # Perform the checks
            pri_check = len(pri_values) >= 3
            sec_check = len(sec_values) >= 1
            
            
            # Additional checks for [Result] section values
            depth_check = depth_value is not None and depth_value >= 1
            rms_check = rms_value is not None and rms_value < 1
            seconds_check = seconds_value is not None and seconds_value < 1
            resp_check = len(resp_values) == 0
            ress_check = len(ress_values) == 0
            mag_check = mag_value is not None and mag_value >= 1
            
            
            return pri_check, sec_check, sta_check, depth_check, rms_check, seconds_check, mag_check, resp_check, ress_check, sta_duplicate, list(set(sta_vol_duplicates)), depth_value, rms_value, seconds_value, mag_value, failed_sta_resp, failed_sta_ress, issuedby_value, sta_used, lat_value, long_value


def find_missing_bulletins():
    try:
        year = entry_fields[0].get()
        month = entry_fields[1].get().zfill(2)
        day = entry_fields[2].get().zfill(2)
            
        if not year.isnumeric() or not month.isnumeric() or not day.isnumeric():
            messagebox.showerror("Invalid Input", "Year, month, and day must be numeric values.")
            return

        if int(month) < 1 or int(month) > 12:
            messagebox.showerror("Invalid Month", "Please enter a valid month (01-12).")
            return

        if int(day) < 1 or int(day) > 31:
            messagebox.showerror("Invalid Day", "Please enter a valid day (01-31).")
            return
            
        eqp_path = f"C:/EARTHQUAKE/SOLUTION_EQ/Sol_{year}/{month}/{day}"
        bul_path = f"C:/EARTHQUAKE/BULLETIN_EQ/Bul_{year}/{month}/{day}"

        eqp_files = [f[:14] for f in os.listdir(eqp_path) if os.path.isfile(os.path.join(eqp_path, f))]
        # Filter files with .html and .jpg extensions
        bul_files = [f[:14] for f in os.listdir(bul_path) if os.path.isfile(os.path.join(bul_path, f)) and (f.endswith('.html'))]

        missing_bulletins = set(eqp_files) - set(bul_files)
        missing_eqp = set(bul_files) - set(eqp_files)

        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
            
        x = glob(f"{eqp_path}/**.eqp")
        
        checks_failed = False
            
        output_text.insert(tk.END, "EQP Checks:\n", "bold")
            
        check_names = ['Primary', 'Secondary', 'Duplicate station(s)', 'Depth', 'RMS', 'SEC', 
                        'Magnitude', 'RESp', 'RESs', 'Station order', 'Existing station']
        

        checks_failed = False
            
        file_number = 0
        for filename in x:
            pri_check, sec_check, sta_check, depth_check, rms_check, seconds_check, mag_check, resp_check, ress_check, list_stations_duplicate, list_vol_stations_duplicates, depth_value, rms_value, seconds_value, mag_value, list_failed_sta_resp, list_failed_sta_ress, name_issuedby, list_stations_used, lat_val, long_val = check_eqp_file(filename)
            failed_checks = []
            tooltip_failed = []
            tooltip_var_dict = {}
            file_name = filename.split("\\")[-1][:14] + filename.split("\\")[-1][22:25]
            # Get distances of station
            dist = st_dist(df, lat_val, long_val)
            
            dst = pd.DataFrame(dist, columns=['Distance'])

            dist_stations = pd.concat([df['Station'].str[:3], dst], axis=1)
            dist_stations.sort_values(by=['Distance'], inplace=True)

            dist_stations.reset_index(drop=True, inplace=True)
                   
            # Dictionary of stations and their distances to the epicenter
            dict_of_sta_and_distances = dist_stations.set_index('Station').T.to_dict('list')
            
            # Convert dict_keys object to a list
            stations_related_to_eq = list(dict_of_sta_and_distances.keys())
            
            stations_func = compare_stations_used_to_quake(list_stations_used, stations_related_to_eq, file_name)
            stations_comparison = stations_func[0]
            result_of_comparison = len(stations_comparison) == 0
            
            list_of_missing_stations = stations_func[1]
            result_of_missing_stations = len(list_of_missing_stations) == 0
            
      
            checks = [pri_check, sec_check, sta_check, depth_check, rms_check, seconds_check, mag_check, resp_check, ress_check, result_of_comparison, result_of_missing_stations]
                
            if not all(checks):
                checks_failed = True
                # Display station and distances from dict_of_sta_and_distances based on stations_comparison
                stations_and_distances = [(station, dict_of_sta_and_distances[station]) for station in stations_comparison if station in dict_of_sta_and_distances]
                
                    
                tooltip_text_list = [
                        "Note: At least three asterisks should remain on the primary side.",
                        "Note: At least one should remain uncrossed on the secondary side.",
                        f"Multiple stations found:\n{', '.join(list_stations_duplicate + list_vol_stations_duplicates)}\n\nNote: For Volcanic stations, at most two similar stations are only accepted.",
                        f"Depth: {depth_value}\nNote: Depth must be greater than or equal to 1",
                        f"RMS: {rms_value}\nNote: RMS must be less than 1",
                        f"SEC: {seconds_value}\nNote: SEC must be less than 1",
                        f"Magnitude: {mag_value}\nNote: Magnitude must be greater than or equal to 1",
                        f"Test failed for station(s):\n{', '.join(list_failed_sta_resp)}\n\nNote: The absolute values of RESp must be less than 1.",
                        f"Test failed for station(s):\n{', '.join(list_failed_sta_ress)}\n\nNote: The absolute values of RESs must be less than 1.",
                        f"Incorrect use of station(s):\n{', '.join(f'{station} ({distance[0]:.2f} km)' for station, distance in stations_and_distances)}\n\nNote: Check the order of stations and their corresponding P-picks.",
                        f"Station(s) not in the database:\n{', '.join(f'{station} in {fn}' for station, fn in list_of_missing_stations)}\n\nNote: Please update the excel file of stations."
                ]
                    
                for i, check in enumerate(checks):
                    if not check:
                            
                        # Append the name of the failed check to the failed_checks list
                        failed_checks.append(check_names[i])
                            
                        tooltip_text = tooltip_text_list[i]
                        tooltip_failed.append(tooltip_text)
   
                output_text.insert(tk.END, f"{file_name} failed in ")
                    

                for i, fail in enumerate(failed_checks):
                    tooltip_var_dict[f"{fail}_{file_number}"] = Tooltip(output_text, tooltip_failed[i])
                    
                    
                for fail, _ in tooltip_var_dict.items():
                        
                    output_text.tag_bind(f"{fail}_test", "<Enter>", lambda event, tooltip=tooltip_var_dict[f"{fail}"]: tooltip.show_tooltip(event))
                    output_text.tag_bind(f"{fail}_test", "<Leave>", lambda event, tooltip=tooltip_var_dict[f"{fail}"]: tooltip.hide_tooltip(event))
                        
                    # Configure the tags to bind tooltips
                    output_text.tag_configure(f"{fail}_test", underline=True)
                        
                    fail_text = fail.split('_')[0]
                    output_text.insert(tk.END, f"{fail_text}", (f"{fail}_test",))
                        
                    if fail_text == failed_checks[-1]:
                        output_text.insert(tk.END, " ")
                    else:
                        output_text.insert(tk.END, ", ")
                        
                if len(failed_checks) > 1:
                    output_text.insert(tk.END, "tests.\n")
                else:
                     output_text.insert(tk.END, "test.\n")
                    
                    
            file_number += 1
                            
                    
        if not checks_failed:
            output_text.insert(tk.END, "All checks passed!\n")
            
        
        output_text.insert(tk.END, "\nMissing Bulletins:\n", "bold")
            
        if missing_bulletins:
                
            output_text.insert(tk.END, "\n".join(missing_bulletins))
        else:
            output_text.insert(tk.END, "No missing bulletins.")
              
        output_text.insert(tk.END, "\n\nMissing EQP:\n", "bold")
            
        if missing_eqp:
            output_text.insert(tk.END, "\n".join(missing_eqp))
        else:
            output_text.insert(tk.END, "No missing EQP.")
                
        output_text.config(state=tk.DISABLED)

    except Exception as e:
            
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"An error occurred: {str(e)}")
        output_text.config(state=tk.DISABLED)
     
    
def extract_eqp_report():
    year = entry_fields[0].get()
    month = entry_fields[1].get().zfill(2)
    day = entry_fields[2].get().zfill(2)
        
    if not year.isnumeric() or not month.isnumeric() or not day.isnumeric():
        messagebox.showerror("Invalid Input", "Year, month, and day must be numeric values.")
        return

    if int(month) < 1 or int(month) > 12:
        messagebox.showerror("Invalid Month", "Please enter a valid month (01-12).")
        return

    if int(day) < 1 or int(day) > 31:
        messagebox.showerror("Invalid Day", "Please enter a valid day (01-31).")
        return
        

    base_path = f"C:/EARTHQUAKE/SOLUTION_EQ/Sol_{year}"
    pathToEvents = os.path.join(base_path, month, day)

    x = glob(f"{pathToEvents}/**.eqp")

    if len(x) == 0:
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, "Folder not found. Check the entered date.")
        output_text.config(state=tk.DISABLED)
    
    else:
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", "end")  # Clear previous output
                
        for filename in x:
            with open(filename) as f:
                lines = f.readlines()

            start_index = next(i for i, line in enumerate(lines) if "[Data]" in line)
            end_index = next(i for i, line in enumerate(lines) if "[Result]" in line)
            mag_start_index = next(i for i, line in enumerate(lines) if "[Magnitude]" in line)
            mag_end_index = next(i for i, line in enumerate(lines) if "[Location]" in line)

            picks = lines[start_index + 2:end_index - 1]
            mags = lines[mag_start_index + 1:mag_end_index - 1]

            count_p = sum(1 for pick in picks if pick[6:13].replace(" ", "") != "")
            count_s = sum(1 for pick in picks if pick[15:22].replace(" ", "") != "")

            eqp_fn = os.path.basename(filename).split("_")
            bul_no = eqp_fn[4].rstrip(".eqp")

                    
            if eqp_fn[4] != eqp_fn[-1]:
                bul = eqp_fn[-1].replace(".eqp", "")
                output_text.insert("end", f"{eqp_fn[2]}Z\t{count_p}\t{count_s}\t{mags[0][3:6]}\t{bul_no}_{bul}\n")
            else:
                output_text.insert("end", f"{eqp_fn[2]}Z\t{count_p}\t{count_s}\t{mags[0][3:6]}\t{bul_no}\n")
                        
        output_text.config(state=tk.DISABLED)
    

def process_dqr():
    try:
        year = entry_fields[0].get()
        month = entry_fields[1].get().zfill(2)
        day = entry_fields[2].get().zfill(2)
        
        base_path = "C:\DQR\STSS_DQR"
        dqr_folder = f"{year}{month}{day}"
        dqr_file = f"DQR_STSS_{year}_{month}_{day}.xlsx"
        dqr_filename = os.path.join(base_path, dqr_folder, dqr_file)
            
        df = pd.read_excel(dqr_filename)
        df.columns = df.iloc[1]
        df.columns.name = None
        df = df[2:]
            
        # Get the directory of the script
        script_directory = os.path.dirname(os.path.abspath(__file__))
        station_excel_path = os.path.join(script_directory, "Stations", 'station_names.xlsx')
            
        station_names = pd.read_excel(station_excel_path)
        station_names['# OF P PICKS'] = station_names['Station'].str[:3].map(df.groupby(df['Station'].str[:3])['P arrival'].count()).fillna(0).astype(int)
        station_names['# OF S PICKS'] = station_names['Station'].str[:3].map(df.groupby(df['Station'].str[:3])['S arrival'].count()).fillna(0).astype(int)
            
        p_sum = station_names['# OF P PICKS'].sum()
        s_sum = station_names['# OF S PICKS'].sum()
            
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"Total P Picks: {p_sum}\nTotal S Picks: {s_sum}", "dqrResult")
        
        output_text.config(state=tk.DISABLED)
            

    except Exception as e:
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"An error occurred: {str(e)}")
        output_text.config(state=tk.DISABLED)


def site_bul():
    year = entry_fields[0].get()
    month = entry_fields[1].get().zfill(2)
    day = entry_fields[2].get().zfill(2)
    time = entry_fields[3].get()

    month_names = {
            '01': 'January',
            '02': 'February',
            '03': 'March',
            '04': 'April',
            '05': 'May',
            '06': 'June',
            '07': 'July',
            '08': 'August',
            '09': 'September',
            '10': 'October',
            '11': 'November',
            '12': 'December'
    }
        
    if not year.isnumeric() or not month.isnumeric() or not day.isnumeric() or not time.isnumeric():
        messagebox.showerror("Invalid Input", "Year, month, day, and time must be numeric values.")
        return

    if month not in month_names:
        messagebox.showerror("Invalid Month", "Please enter a valid month (01-12).")
        return

    if int(day) < 1 or int(day) > 31:
        messagebox.showerror("Invalid Day", "Please enter a valid day (01-31).")
        return
        
    output_text.config(state=tk.NORMAL)
    output_text.delete("1.0", tk.END)

    result = f"For {time}Z:\n"
    found_bul = False
    base_site = f"https://earthquake.phivolcs.dost.gov.ph/{year}_Earthquake_Information/{month_names[month]}/{year}_{month}{day}_"
    bul_html = ['B1', 'B1F', 'B2', 'B2F', 'B3', 'B3F', 'B4', 'B4F', 'B5', 'B5F', 'B6', 'B6F']
    for html in bul_html:
        url = f"{base_site}{time}_{html}.html"
        response = rq.get(url, verify=False)
        if response.ok:
            result += f"Has {html}\n"
            found_bul = True
        
    if not found_bul:
        result = 'No bulletin'
        
    output_text.insert(tk.END, result, "dqrResult")
    output_text.config(state=tk.DISABLED)
        

def run_selenium_script(year, month, day):
    global driver  # Access the global driver instance

    try:
        if driver is None:
            # Chrome driver setup
            driver_path = os.path.join(script_directory, "driver", "chromedriver.exe")
            ser = Service(driver_path)
            op = webdriver.ChromeOptions()
            op.add_argument("--start-maximized")
            driver = webdriver.Chrome(service=ser, options=op)

        # Your Selenium script here
        date = datetime(year, month, day).strftime('%Y-%m-%d')
        driver.get("https://earthquake.usgs.gov/earthquakes/search/")
        start_time = driver.find_element(By.NAME, "starttime")
        start_time.clear()
        start_time.send_keys(date + " 00:00:00")

        end_time = driver.find_element(By.NAME, "endtime")
        end_time.clear()
        end_time.send_keys(date + " 23:59:59")

        toggle = driver.find_element(By.CLASS_NAME, "toggle")
        driver.execute_script("arguments[0].setAttribute('class','toggle toggle-visible')", toggle)

        driver.find_element(By.NAME, "maxlatitude").send_keys("54.47004")
        driver.find_element(By.NAME, "minlatitude").send_keys("-39.63954")
        driver.find_element(By.NAME, "maxlongitude").send_keys("-111.26953")
        driver.find_element(By.NAME, "minlongitude").send_keys("-308.14453")
        driver.find_element(By.ID, "fdsn-submit").click()

        # Open geofon search in a new window
        driver.execute_script("window.open('https://geofon.gfz-potsdam.de/eqinfo/form.php?lang=en')")
        driver.switch_to.window(driver.window_handles[-1])

        # Set up geofon
        driver.find_element(By.NAME, "datemin").send_keys(date)
        driver.find_element(By.NAME, "datemax").send_keys(date)
        driver.find_element(By.NAME, "latmax").send_keys("54.47004")
        driver.find_element(By.NAME, "latmin").send_keys("-39.63954")
        driver.find_element(By.NAME, "lonmax").send_keys("-111.26953")
        driver.find_element(By.NAME, "lonmin").send_keys("-308.14453")
        driver.find_element(By.NAME, "nmax").send_keys("500")
        driver.find_element(By.CLASS_NAME, "btn-primary").click()
            
    except Exception as e:
        messagebox.showerror("Error", f"An error occurred: {str(e)}")


def execute_script():
    year_input = entry_fields[0].get()
    month_input = entry_fields[1].get()
    day_input = entry_fields[2].get()

    if year_input and month_input and day_input:
        try:
            year = int(year_input)
            month = int(month_input)
            day = int(day_input)
            run_selenium_script(year, month, day)
        except ValueError:
            messagebox.showerror("Error", "Invalid input. Year, Month, and Day must be integers.")
            return
    else:
        messagebox.showerror("Error", "Please fill in all fields.")


def gsheets_tool():
    def get_last_row(worksheet):
        return f"{len(list(filter(None, worksheet.col_values(1)))) + 1}"

    def count_digits(s):
        return sum(c.isdigit() for c in s)

    def show_popup():
        # Display a popup message when the processing is done
        messagebox.showinfo("Processing Complete", "The data processing is complete!")

    def get_other_sheet_names(current_sheet_name, num_days_before, num_days_after):
        current_date = datetime.strptime(current_sheet_name, "%Y%m%d")
        other_sheet_names = []
        for i in range(num_days_before, 0, -1):
            previous_date = current_date - timedelta(days=i)
            previous_sheet_name = previous_date.strftime("%b%d")  #  Sheet name date format
            other_sheet_names.append(previous_sheet_name)
        other_sheet_names.append(current_sheet_name)
        for i in range(1, num_days_after + 1):
            next_date = current_date + timedelta(days=i)
            next_sheet_name = next_date.strftime("%b%d")   #  Sheet name date format
            other_sheet_names.append(next_sheet_name)
        return other_sheet_names

    def get_values(worksheet, sheet_name):
        try:
           workbook = load_workbook(worksheet)
           sheet = workbook[sheet_name]
        except KeyError:
           return None, None, None
       
        column = sheet['A']
        last_row = 0
        consecutive_empty = 0
        for cell in column:
            if cell.value is not None:
                consecutive_empty = 0
                last_row = cell.row
            else:
                consecutive_empty += 1
                if consecutive_empty > 1:
                    break
        
        column_values = []
        range_str = f"D4:D{last_row}"
        for index, row in enumerate(sheet[range_str], start=4):
            for cell in row:
                column_values.append((index, cell.value))
                
        return sheet, column_values, last_row

    def process_spreadsheet(link, sheet_name, num_days_before, num_days_after, progress_queue, lock):
        global progress_window
        try:
            # Attempt to connect to the spreadsheet URL
            while True:
                try:
                    # Calculate sheet_date from the input sheet_name
                    sheet_date = datetime.strptime(sheet_name, "%Y%m%d")
                    current_sheet_name = sheet_date.strftime("%Y%m%d")
                    g_account = gspread.service_account(credentials_path)
                    gsheets1 = g_account.open_by_url(link)
                    break  # Break the loop if connection successful
                except Exception as e:
                    messagebox.showinfo("Error", f"Connection error: {str(e)}")
                    
                    time.sleep(5)  # Wait for 5 seconds before the next attempt

            # Connect to the main worksheet of the first spreadsheet
            main_worksheet = gsheets1.worksheet(current_sheet_name)
            main_worksheet.delete_rows(2)
            main_worksheet.sort((4, 'asc'), range=f"A3:J{get_last_row(main_worksheet)}")

            data = main_worksheet.get_all_values()

            year = sheet_date.strftime("%Y")
            month = sheet_date.strftime("%m")

            worksheet_file = f"{year}_{month}_phase data.xlsx"
            worksheet_path = os.path.join(phase_data_folder, worksheet_file)


            new_data = []
            prev_p_time = None
            stations = []
            data_error = []

            for row in data[2:]:
                if count_digits(row[3][:6]) < 6:
                    data_error.append(row)
                    continue
                try:
                    p_time = datetime.strptime(row[3], '%H%M%S.%f')
                except ValueError:
                    data_error.append(row)
                    continue
                sta = row[0]
                if prev_p_time is None or (p_time - prev_p_time).total_seconds() > 40 or sta in stations:
                    new_data.append([''] * len(data[0]))
                    stations = []
                new_data.append(row)
                prev_p_time = p_time
                stations.append(sta)

            main_worksheet.clear()
            main_worksheet.append_rows(data[:2])
            set_frozen(main_worksheet, rows=2)

            # get_other_sheet_names(current_sheet_name, num_days_before, num_days_after):
            other_sheet_names = get_other_sheet_names(current_sheet_name, num_days_before, num_days_after)

            
            total_rows = len(new_data) + len(data_error)
            progress_step = 100 / total_rows
            current_progress = 0

            
            for ind, row in enumerate(data_error, start=2):
                main_worksheet.append_row(row, table_range=f"A{ind + 1}")
                main_worksheet.format(f"D{ind + 1}", {"backgroundColor": {"red": 1, "green": 0.51, "blue": 0.5}})
                time.sleep(1)
                current_progress += progress_step
                with lock:
                    progress_queue.put(current_progress)
            
            ind = len(data_error) + 2

            for row in new_data:
                p_arrival = row[3]  # P arrival in the 4th column
                sta = row[0]
                if sta == '':
                    main_worksheet.insert_row([], ind + 1)
                else:
                    duplicate_found = False
                    found_in = []
                    for sheet_names in other_sheet_names:
                        if sheet_names != current_sheet_name:
                            sheets, column_values, _ = get_values(worksheet_path, sheet_names)
                            if column_values is None:  # Sheet not found, continue to the next iteration
                                continue
                            for index, value in column_values:
                                first_column_value = sheets[f"A{index}"].value
                                if p_arrival == value and sta == first_column_value:
                                    duplicate_found = True
                                    found_in.append(sheet_names)
                                    
                    if duplicate_found:
                        main_worksheet.append_row(row, table_range=f"A{ind + 1}")
                        found_in_list = ", ".join(found_in)
                        main_worksheet.update(f"J{ind + 1}", f"Found in: {found_in_list}")
                        main_worksheet.format(f"J{ind + 1}", {"backgroundColor": {"red": 1, "green": 0.51, "blue": 0.5}})
                    else:
                        main_worksheet.append_row(row, table_range=f"A{ind + 1}")
                    
                ind += 1
                
                # Update the progress queue with the current progress
                current_progress += progress_step
                with lock:
                    progress_queue.put(current_progress)
                
                time.sleep(1)
            
            val = main_worksheet.acell('A3').value
            
            if not val:
                main_worksheet.delete_rows(3)
            
            check_gsheets_last_row = int(get_last_row(main_worksheet))
            
            current_dqr = datetime.strptime(current_sheet_name, "%Y%m%d")
            dqr_name = current_dqr.strftime("%b%d")
            *_, dqr_last_row = get_values(worksheet_path, dqr_name)
            dqr_last_row = int(dqr_last_row)
            
          
            if check_gsheets_last_row == dqr_last_row:

                autogsheets.after(100, show_popup)  # Schedule the popup message after a short delay
            
                # Close the progress window
                progress_window.destroy()
    
                # After processing is complete, restore the main window
                autogsheets.deiconify()
                
            else:
               messagebox.showinfo("Error", "DQR and Google Sheets rows are not equal. Please repeat the process.")
               # Close the progress window
               progress_window.destroy()
               
               # After processing is complete, restore the main window
               autogsheets.deiconify()
            
        except Exception as e:
            messagebox.showinfo("Error", f"An error occured: {str(e)}. Please repeat the process")
            # Close the progress window
            progress_window.destroy()
            
            # After processing is complete, restore the main window
            autogsheets.deiconify()
            
    # Function to update the progress bar in the main thread
    def update_progress(progress_queue, lock):
        global progress_var, progress_window
        try:
            while True:
                with lock:
                    value = progress_queue.get_nowait()
                    progress_var.set(value)
                    progress_window.update_idletasks()
        except queue.Empty:
            autogsheets.after(100, update_progress, progress_queue, lock )
            
    def start_processing():
        global progress_window, progress_var
        link = link_entry.get()
        sheet_name = sheet_name_entry.get()
        num_days_before = 5
        num_days_after = 5

        # Disable the Start Processing button while processing
        start_button.config(state=tk.DISABLED)

        # Hide the main window while processing
        autogsheets.withdraw()

        # Create a separate window for the progress bar
        progress_window = tk.Toplevel(autogsheets)
        progress_window.iconbitmap(icon_path)
        progress_window.title("Processing...")
        progress_window.geometry("300x75")

        # Center the progress window on the screen
        center_window(progress_window)

        # Create a progress bar in the progress window
        progress_var = tk.DoubleVar()
        progress_bar = ttk.Progressbar(progress_window, mode='determinate', variable=progress_var, maximum=100)
        progress_bar.pack(pady=20)
        
        # Create a lock to protect the progress queue
        lock = threading.Lock()

         # Create a separate thread for processing
        processing_thread = threading.Thread(target=process_spreadsheet, args=(link, sheet_name, num_days_before, num_days_after, progress_queue, lock))
        
        # Start the progress bar in the main thread
        autogsheets.after(100, update_progress, progress_queue, lock)

        # Check the progress in a separate thread
        autogsheets.after(100, check_progress, processing_thread)
        
        processing_thread.start()


    def check_progress(processing_thread):
        
        if processing_thread.is_alive():
            # Continue checking progress after a short delay
            autogsheets.after(100, check_progress, processing_thread)
        else:
            start_button.config(state=tk.NORMAL)
            
    progress_queue = queue.Queue()
    

    # Construct absolute paths based on the script directory
    icon_path = os.path.join(script_directory, "icons", "gsheets_icon.ico")
    gsheet_path = os.path.join(script_directory, "gsheet_credentials")
    # List files in the folder
    folder_contents = os.listdir(gsheet_path)
    
    # Filter for JSON files in the folder
    credential_files = [f for f in folder_contents if f.endswith(".json")]
    json_file_name = credential_files[0]
    credentials_path = os.path.join(gsheet_path, json_file_name)
    
    phase_data_folder = os.path.join(script_directory, "PhaseData")
    
    # Create the main application window
    autogsheets = tk.Toplevel(root)
    autogsheets.protocol("WM_DELETE_WINDOW", lambda ags=autogsheets: on_close(ags))
    autogsheets.title("Google Sheets Processing")
    autogsheets.geometry("440x110")
    
    # Hide the main root window when the Toplevel is opened
    hide_root()
    
    # Center the main window on the screen
    center_window(autogsheets)
    autogsheets.minsize(440, 110)
    autogsheets.maxsize(440, 110)
    autogsheets.iconbitmap(icon_path)

    # Labels
    link_label = tk.Label(autogsheets, text="Google Sheets Link:")
    link_label.grid(row=0, column=0, padx=5, pady=5, sticky=tk.W)
    sheet_name_label = tk.Label(autogsheets, text="Sheet Name (YYYYMMDD):")
    sheet_name_label.grid(row=1, column=0, padx=5, pady=5, sticky=tk.W)


    # Entry fields with validation
    link_entry = tk.Entry(autogsheets, width=40)
    link_entry.grid(row=0, column=1, padx=5, pady=5)
    sheet_name_entry = tk.Entry(autogsheets, width=40)
    sheet_name_entry.grid(row=1, column=1, padx=5, pady=5)


    # Button
    start_button = tk.Button(autogsheets, text="Start Processing", command=start_processing)
    start_button.grid(row=4, column=0, columnspan=2, padx=5, pady=10)


def display_about_dialog():
    about_text = (
        "STDM Tools\n"
        "Version 2.2\n"
        "Created by: ABD\n"
        "Description: This application is designed to assist with earthquake analysis.\n"
        
    )
    messagebox.showinfo("About", about_text)
    

# Function to open a tool using subprocess
def open_tool(tool_function):
    try:
        tool_function()
    except Exception as e:
        messagebox.showerror("Error", f"Error running tool: {str(e)}")

# Function to hide the main root window
def hide_root():
    global root_hidden
    root.withdraw()
    root_hidden = True

# Function to restore the main root window
def restore_root():
    global root_hidden
    root.deiconify()
    root_hidden = False
    
def on_close(window):
    # Exit 
    if isinstance(window, tk.Tk):
        window.destroy()
        sys.exit()
    else:
        window.destroy()
        restore_root()

def center_window(window):
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    x_offset = (window.winfo_screenwidth() - width) // 2
    y_offset = (window.winfo_screenheight() - height) // 2
    window.geometry(f"+{x_offset}+{y_offset}")

def update_resolution(event):
    selected_tab = event.widget.select()
    tab_index = event.widget.index(selected_tab)
    
    if tab_index == 0:
        # Resolution settings for Tab 1
        root.geometry("200x100")
        root.maxsize(200, 100)
        root.minsize(200, 100)
    elif tab_index == 1:
        # Resolution settings for Tab 2
        root.geometry("470x385")
        root.maxsize(470, 385)
        root.minsize(470, 385)

    center_window(root)

# Function to ensure only 1 instance of application is running
def is_another_instance_running():
    mutex_name = "Global\\MyAppMutex"  # Unique mutex name
    mutex = CreateMutex(None, True, mutex_name)
    if GetLastError() == 183:  # ERROR_ALREADY_EXISTS
        CloseHandle(mutex)
        return True
    else:
        return False

if __name__ == '__main__':
     
    # Initialize a global driver instance
    driver = None
    processing_thread = None
    root_hidden = False
    
    if is_another_instance_running():
        MessageBox(None, "Another instance is already running.", "Error", 0 | 0x30)
        sys.exit(1)
        
    else:
        from tkinter import Scrollbar, messagebox, ttk, Menu
        import tkinter as tk
        
        # Create the main menu window
        root = tk.Tk()
        root.title("STDM Tools")
       
        # Get the directory of the script
        if getattr(sys, 'frozen', False):
        # Running as a bundled executable
            script_directory = sys._MEIPASS
        else:
            # Running as a script
            script_directory = os.path.abspath(os.path.dirname(__file__))
            
        station_dir = os.path.join(script_directory, 'Stations', 'station_names.xlsx')
        df = pd.read_excel(station_dir)
        
        # Construct absolute paths based on the script directory
        icon_path = os.path.join(script_directory, "icons", "tools_icon.ico")
        root.iconbitmap(icon_path)
        
        main_menu = Menu(root)
        root.config(menu=main_menu)
        
        # Add the "About" option to the main menu bar
        main_menu.add_command(label="About", command=display_about_dialog)
        
        # Create a tabbed interface
        tab_control = ttk.Notebook(root)
    
        # Create tabs for Main Tools and Report Tools
        main_tools_tab = ttk.Frame(tab_control)
        check_tools_tab = ttk.Frame(tab_control)
    
        # Add the tabs to the tabbed interface
        tab_control.add(main_tools_tab, text="Main Tool")
        tab_control.add(check_tools_tab, text="Other Tools")
    
        # Define the tools for each category
        main_tools =  {"name": "Automate Google Sheets", "function": gsheets_tool}
        
        
        check_report_tools = [
            {"name": "Search USGS and GEOFON", "function": execute_script},
            {"name": "Search website bulletin", "function": site_bul},
            {"name": "Check EQP & Bulletins", "function": find_missing_bulletins},
            {"name": "Extract EQP", "function": extract_eqp_report},
            {"name": "Count P and S", "function": process_dqr}
        ]
    
        # Create entry fields
        entry_fields = []
    
        # Date labels
        labels = ["Year", "Month", "Day", "Time"]
    
        for i in range(len(labels)):
            # Create label
            label = ttk.Label(check_tools_tab, text=f"{labels[i]}: ")
            label.grid(row=i+1, column=0, pady=5, sticky=tk.E)
    
            # Create entry field
            entry = ttk.Entry(check_tools_tab)
            entry.grid(row=i+1, column=1, pady=5, sticky=tk.W)
            entry_fields.append(entry)
            entry.config(validate="key", validatecommand=(check_tools_tab.register(lambda x: x.isnumeric() or x == ""), "%P"))
    
        # Main tool button
        maintool_button = tk.Button(main_tools_tab, text=main_tools["name"], command=lambda func=main_tools["function"]: open_tool(func))
        maintool_button.pack(pady=20)
    
        # Create buttons for each tool in the Main Tools category
        for i, tool in enumerate(check_report_tools):
    
            tool_button = tk.Button(check_tools_tab, text=tool["name"], command=lambda func=tool["function"]: open_tool(func))
            tool_button.grid(row=i, column=2, pady=5, sticky=tk.W)
    
        # Bind tab change event to update_resolution function
        tab_control.bind("<<NotebookTabChanged>>", update_resolution)
        
        # Pack the tabbed interface
        tab_control.pack()
    
        # Create a Text widget for output in Tab 2
        output_text = tk.Text(check_tools_tab, height=10, width=50, state=tk.DISABLED)
        output_text.grid(row=5, columnspan=3, padx=5, pady=5)
        output_text.tag_configure("bold", font=("Helvetica", 11, "bold"))
        output_text.tag_configure("dqrResult", font=("Helvetica", 15))
    
        # Scrollbar for the output_text
        scrollbar = Scrollbar(check_tools_tab, command=output_text.yview)
        scrollbar.grid(row=5, column=3, sticky="ns")
        output_text.config(yscrollcommand=scrollbar.set)
    
        # Specify what should happen when the window is closed
        root.protocol("WM_DELETE_WINDOW", lambda rt=root: on_close(rt))
        root.mainloop()

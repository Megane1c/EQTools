import os
import sys
import ctypes
import pandas as pd
from numpy import arcsin, sin, cos, sqrt
from math import pi
from glob import glob
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from threading import Timer

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

# Observer with debouncer
class FileModifiedHandler(FileSystemEventHandler):
    
    def __init__(self, debounce_delay):
        super().__init__()
        self.debounce_delay = debounce_delay
        self.last_event = None
        self.debouncer = None

    def process_file(self):
        if self.last_event is not None:
            # Process the file
            process_eqp_files(self.last_event.src_path)
    
    def debounce(self):
        if self.debouncer:
            self.debouncer.cancel()
        self.debouncer = Timer(self.debounce_delay, self.process_file)
        self.debouncer.start()
    
    def on_created(self, event):
        self.last_event = event
        self.debounce()
        
    def on_modified(self, event):
        self.last_event = event
        self.debounce()
        

# Calculates the distance of a station to an epicenter
def st_dist(df, x, y):
    dist = 6371*2*arcsin(sqrt(sin((df['Latitude']*pi/180 - x*pi/180)/2)**2 + 
                        cos(x*pi/180)*cos(df['Latitude']*pi/180)*sin((df['Longitude']*pi/180 - y*pi/180)/2)**2))
    
    return dist
            

# For selecting mode
def toggle_input_method():
    global selected_method
    selected_method = input_method.get()
    
    missing_bul.grid_rowconfigure(7, weight=1)
    
    if selected_method == 1:
        # Individual date input method
        year_label.grid()
        year_entry.grid()
        month_label.grid()
        month_entry.grid()
        day_label.grid()
        day_entry.grid()
        find_button.grid()
        select_folders_button.grid_remove()  # Hide the folder selection button
        startprocess_button.grid_remove() # Hide the startprocess button
        endprocess_button.grid_remove() # Hide the endprocess button
    
    elif selected_method == 2 or selected_method == 3:
        # Folder selection input method
        year_label.grid_remove()  # Hide year label
        year_entry.grid_remove()  # Hide year entry
        month_label.grid_remove()  # Hide month label
        month_entry.grid_remove()  # Hide month entry
        day_label.grid_remove()  # Hide day label
        day_entry.grid_remove()  # Hide day entry
        
        find_button.grid_remove()
        startprocess_button.grid_remove()
        endprocess_button.grid_remove()
        
        select_folders_button.grid()  # Show the folder selection button
        
    elif selected_method == 4:
        # Individual date input method
        missing_bul.grid_rowconfigure(7, weight=0)
        year_label.grid()
        year_entry.grid()
        month_label.grid()
        month_entry.grid()
        day_label.grid()
        day_entry.grid()
        startprocess_button.grid()
        
        find_button.grid_remove() # Hide the find button
        select_folders_button.grid_remove()  # Hide the folder selection button
        

    if output_text.cget("state") == "disabled":
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
   
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


def start_process():
    global observer, selected_method
    # Determine the folder to watch based on the selected radio button
    year = year_entry.get()
    month = month_entry.get().zfill(2)
    day = day_entry.get().zfill(2)
        
    if not year.isnumeric() or not month.isnumeric() or not day.isnumeric():
        messagebox.showerror("Invalid Input", "Year, month, and day must be numeric values.")
        return

    if int(month) < 1 or int(month) > 12:
        messagebox.showerror("Invalid Month", "Please enter a valid month (01-12).")
        return

    if int(day) < 1 or int(day) > 31:
        messagebox.showerror("Invalid Day", "Please enter a valid day (01-31).")
        return
        
    eqp_path = f"C:\EARTHQUAKE\SOLUTION_EQ\Sol_{year}\{month}\{day}"
    if not os.path.exists(eqp_path):
        os.makedirs(eqp_path)
        
    endprocess_button.grid()
    startprocess_button.grid_remove()
    
    individual_date_radio.config(state=tk.DISABLED)
    folder_selection_radio.config(state=tk.DISABLED)
    generate_report_radio.config(state=tk.DISABLED)
    
    year_entry.config(state=tk.DISABLED)
    month_entry.config(state=tk.DISABLED)
    day_entry.config(state=tk.DISABLED)
    
    # Configure the observer if it's not configured already
    if not observer:
        observer_handler = FileModifiedHandler(1)
        observer = Observer()
        observer.schedule(observer_handler, eqp_path, recursive=True)
        observer.start()
        
def end_process():
    global observer
    # Show the start button
    startprocess_button.grid()
    # Hide the end button
    endprocess_button.grid_remove()
    
    # Enable other buttons
    individual_date_radio.config(state=tk.NORMAL)
    folder_selection_radio.config(state=tk.NORMAL)
    generate_report_radio.config(state=tk.NORMAL)
    
    year_entry.config(state=tk.NORMAL)
    month_entry.config(state=tk.NORMAL)
    day_entry.config(state=tk.NORMAL)
    
    if output_text.cget("state") == "disabled":
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        
    if observer:
        observer.stop()
        observer.join()
        observer = None
        
        
def generate_text_report(data):
    
    
    if len(data) != 0:
        report_dir = os.path.join(script_directory, "Reports")
        # Create the "Reports" folder if it doesn't exist
        if not os.path.exists(report_dir):
            os.makedirs(report_dir)
        
        # Define tab widths based on the header
        issued_by_width = 15
        eqp_file_width = 20
        error_report_width = 50
    
        # Group data by year and month for file naming
        file_names = {}
        for entry in data:
            year_month = entry["eqp_file"].split("_")[0] + "_" + entry["eqp_file"].split("_")[1][:2]
            if year_month not in file_names:
                file_names[year_month] = []
            file_names[year_month].append(entry)
    
        # Write error reports for each year and month group
        for year_month, entries in file_names.items():
            # Generate file name based on year and month
            file_name = f"{year_month}-Report.txt"
            
            file_path = os.path.join(report_dir, file_name)
    
            # Open the text file for writing (or appending if file exists)
            with open(file_path, "a") as file:
                # Write the column headers with specific tab widths
                if file.tell() == 0:  # Check if file is empty
                    file.write("FOR CORRECTIONS\n")
                    file.write(f"Issued By{' '*(issued_by_width-len('Issued By'))}\t"
                               f"EQP File{' '*(eqp_file_width-len('EQP File'))}\t"
                               f"Error Report{' '*(error_report_width-len('Error Report'))}\n")
    
                # Write data for each entry with aligned tabbing
                for entry in entries:
                    issued_by = entry["name"]
                    eqp_file = entry["eqp_file"]
                    error_reports = entry["errors"]
    
                    # Write the first report for each entry
                    if error_reports:
                        first_report = error_reports[0]
                        file.write(f"{issued_by:{issued_by_width}}\t"
                                   f"{eqp_file:{eqp_file_width}}\t"
                                   f"{first_report:{error_report_width}}\n")
                        # Write additional reports below the first one
                        for report in error_reports[1:]:
                            file.write(f"{' '*issued_by_width}\t"
                                       f"{' '*eqp_file_width}\t"
                                       f"{report:{error_report_width}}\n")
                    else:
                        file.write(f"{issued_by:{issued_by_width}}\t"
                                   f"{eqp_file:{eqp_file_width}}\n")
                    
        messagebox.showerror("Report", "Errors were found.")


def process_eqp_files(current_file=None):
    global selected_method
    
    if output_text.cget("state") == "disabled":
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        
    try:
        selected_method = input_method.get()
        if selected_method == 1:
            year = year_entry.get()
            month = month_entry.get().zfill(2)
            day = day_entry.get().zfill(2)
                
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
                
            if os.path.exists(eqp_path):
                x = glob(f"{eqp_path}/**.eqp")
            else:
                messagebox.showerror("Invalid Folder", "Double check the date to be processed.")
                return
                
        elif selected_method == 2 or selected_method == 3:
            # Folder selection input method
            
            eqp_folder_path = filedialog.askdirectory(title='Select a Folder')
            
            if os.path.exists(eqp_folder_path):
            
                # Create an empty list to store the paths of all ".eqp" files
                x = []
                    
                # Walk through the directory tree and find ".eqp" files
                for root, _, files in os.walk(eqp_folder_path):
                    for file in files:
                        if file.endswith(".eqp"):
                            x.append(os.path.join(root, file))
                            
            else:
                messagebox.showerror("Invalid Folder", "Please select a folder.")
                return
                        
        else:
            # For selected method = 3
            x = [current_file]

        check_names = ['Primary', 'Secondary', 'Duplicate station(s)', 'Depth', 'RMS', 'SEC', 
                        'Magnitude', 'RESp', 'RESs', 'Station order', 'Existing station']
        
        data = []

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
            if selected_method == 4:
                output_text.insert(tk.END, "RESULT\n\n", "bold")
            
            if not all(checks):
                checks_failed = True
                # Display station and distances from dict_of_sta_and_distances based on stations_comparison
                stations_and_distances = [(station, dict_of_sta_and_distances[station]) for station in stations_comparison if station in dict_of_sta_and_distances]
                
                if selected_method == 3:
                
                    generate_report_text_list = [
                            "At least three should remain on the primary side.",
                            "At least one should remain uncrossed on the secondary side.",
                            f"Duplicate stations: {', '.join(list_stations_duplicate + list_vol_stations_duplicates)}",
                            f"Depth: {depth_value}",
                            f"RMS: {rms_value}",
                            f"SEC: {seconds_value}",
                            f"Magnitude: {mag_value}",
                            f"RESp failed for station(s): {', '.join(list_failed_sta_resp)}",
                            f"RESs failed for station(s): {', '.join(list_failed_sta_ress)}",
                            f"Incorrect use of station(s): {', '.join(f'{station} ({distance[0]:.2f} km)' for station, distance in stations_and_distances)}",
                            f"Station(s) not in the database: {', '.join(f'{station} in {fn}' for station, fn in list_of_missing_stations)}"
                    ]
                
                else:
                    
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
                        
                        if selected_method == 3:
                            tooltip_text = generate_report_text_list[i]
                        else:
                            tooltip_text = tooltip_text_list[i]
                            
                        tooltip_failed.append(tooltip_text)
                        
             
                if selected_method == 3:
                    data.append({"name": name_issuedby, "eqp_file": file_name, "errors": tooltip_failed})
                
                
                if selected_method == 4:
                    file_name = file_name[10:]
                    
                if selected_method != 3:
                    output_text.insert(tk.END, f"{file_name} failed in ", "normal")
                    
                    for i, fail in enumerate(failed_checks):
                        tooltip_var_dict[f"{fail}_{file_number}"] = Tooltip(output_text, tooltip_failed[i])
                        
                        
                    for fail, _ in tooltip_var_dict.items():
                            
                        output_text.tag_bind(f"{fail}_test", "<Enter>", lambda event, tooltip=tooltip_var_dict[f"{fail}"]: tooltip.show_tooltip(event))
                        output_text.tag_bind(f"{fail}_test", "<Leave>", lambda event, tooltip=tooltip_var_dict[f"{fail}"]: tooltip.hide_tooltip(event))
                            
                        # Configure the tags to bind tooltips
                        output_text.tag_configure(f"{fail}_test", underline=True, font=("Helvetica", 11))
                            
                        fail_text = fail.split('_')[0]
                        output_text.insert(tk.END, f"{fail_text}", (f"{fail}_test",))
                            
                        if fail_text == failed_checks[-1]:
                            output_text.insert(tk.END, " ")
                        else:
                            output_text.insert(tk.END, ", ")
                            
                    if len(failed_checks) > 1:
                        output_text.insert(tk.END, "tests.\n\n", "normal")
                    else:
                        output_text.insert(tk.END, "test.\n\n", "normal")
                    
                    
            file_number += 1
                                  
        if not checks_failed:
            if selected_method == 4:
                name_of_file = x[0].split("\\")[-1][10:14] + x[0].split("\\")[-1][22:25]
                output_text.insert(tk.END, f"{name_of_file}", "normal")
                output_text.insert(tk.END, " passed\n\n", "bold")
                
            else:
                output_text.insert(tk.END, "All files passed!\n", "bold")
            
        
        if selected_method == 3 and checks_failed:
            generate_text_report(data)
            output_text.insert(tk.END, "Please find the ", "normal")
            output_text.insert(tk.END, "Reports ", "bold")
            output_text.insert(tk.END, "folder to see the results.\n", "normal")
            
            

        output_text.config(state=tk.DISABLED)
        
    except Exception as e:
        output_text.config(state=tk.NORMAL)
        output_text.delete("1.0", tk.END)
        output_text.insert(tk.END, f"An error occurred: {str(e)}.")
        output_text.config(state=tk.DISABLED)


# Function to ensure only 1 instance of application is running
def is_another_instance_running():
    mutex_name = "Global\\EQPCheckerAppMutex"  # Unique mutex name
    mutex = CreateMutex(None, True, mutex_name)
    if GetLastError() == 183:  # ERROR_ALREADY_EXISTS
        CloseHandle(mutex)
        return True
    else:
        return False
    

def display_about_dialog():
    about_text = (
        "STDM EQP Checker\n"
        "Version 1.0\n"
        "Created by: ABD\n"        
    )
    messagebox.showinfo("About", about_text)


if __name__ == '__main__':
    
    if is_another_instance_running():
        MessageBox(None, "Another instance is already running.", "Error", 0 | 0x30)
        sys.exit(1)
        
    else:
        import tkinter as tk
        from tkinter import messagebox, filedialog, Text, Scrollbar, Menu
        
        
        # Create the main window
        missing_bul = tk.Tk()
        missing_bul.title("EQP Checker")
        
        # Add the "About" option to the main menu bar
        main_menu = Menu(missing_bul)
        missing_bul.config(menu=main_menu)
        
        main_menu.add_command(label="About", command=display_about_dialog)
        
        
        # Global variable to track the selected input method
        selected_method = tk.IntVar()
        
        # Add radio buttons to select input method
        input_method_label = tk.Label(missing_bul, text="Select Input Method:")
        input_method_label.grid(row=0, column=0, pady=5)
        
        input_method = tk.IntVar()
        input_method.set(1)  # Default to the first method (Individual date input)
        
        individual_date_radio = tk.Radiobutton(missing_bul, text="Date                           ", variable=input_method, value=1, command=toggle_input_method)
        individual_date_radio.grid(row=0, column=1, pady=5)
        
        folder_selection_radio = tk.Radiobutton(missing_bul, text="Folder Selection       ", variable=input_method, value=2, command=toggle_input_method)
        folder_selection_radio.grid(row=1, column=1, pady=5)
        
        generate_report_radio = tk.Radiobutton(missing_bul, text="Generate Report       ", variable=input_method, value=3, command=toggle_input_method)
        generate_report_radio.grid(row=2, column=1, pady=5)
    
        autodetect_radio = tk.Radiobutton(missing_bul, text="Check continuously ", variable=input_method, value=4, command=toggle_input_method)
        autodetect_radio.grid(row=3, column=1, pady=5)
        
        # Labels and Entry fields for individual date input
        year_label = tk.Label(missing_bul, text="Year:")
        year_entry = tk.Entry(missing_bul)
        month_label = tk.Label(missing_bul, text="Month:")
        month_entry = tk.Entry(missing_bul)
        day_label = tk.Label(missing_bul, text="Day:")
        day_entry = tk.Entry(missing_bul)
        
        # Arrange the widgets on the grid
        year_label.grid(row=4, column=0, pady=5)
        year_entry.grid(row=4, column=1, pady=5)
        month_label.grid(row=5, column=0, pady=5)
        month_entry.grid(row=5, column=1, pady=5)
        day_label.grid(row=6, column=0, pady=5)
        day_entry.grid(row=6, column=1, pady=5)
        
        ## BUTTONS ##
        # Button to trigger folder selection
        select_folders_button = tk.Button(missing_bul, text="Select a Folder", command=process_eqp_files)
        select_folders_button.grid(row=6, column=0, columnspan=2, pady=5)
        select_folders_button.grid_remove()  # Initially hide the folder selection button
        
        # Button to trigger finding missing bulletins
        find_button = tk.Button(missing_bul, text="Process", command=process_eqp_files)
        find_button.grid(row=7, columnspan=2, pady=10)
    
        # Buttons for checking continuously
        startprocess_button = tk.Button(missing_bul, text="START", command=start_process)
        startprocess_button.grid(row=7, columnspan=2, pady=10)
        startprocess_button.grid_remove()  # Initially hide the startprocess button
    
        endprocess_button = tk.Button(missing_bul, text="END", command=end_process)
        endprocess_button.grid(row=7, columnspan=2, pady=10)
        endprocess_button.grid_remove()  # Initially hide the endprocess button
    
        # Configure row 7 to expand when the window is resized
        missing_bul.grid_rowconfigure(8, weight=1)
        
        # Output Text
        output_text = Text(missing_bul, height=10, width=50)
        output_text.grid(row=8, columnspan=2, padx=10, pady=5, sticky="nsew")
    
        # Scrollbar for the output_text
        scrollbar = Scrollbar(missing_bul, command=output_text.yview)
        scrollbar.grid(row=8, column=2, sticky="ns")
    
        output_text.config(yscrollcommand=scrollbar.set)
        output_text.config(state=tk.DISABLED)
        output_text.tag_configure("bold", font=("Helvetica", 11, "bold"))
        output_text.tag_configure("normal", font=("Helvetica", 11))
        
        width = 450
        height = 390
    
        # Get the screen width and height
        screen_width = missing_bul.winfo_screenwidth()
        screen_height = missing_bul.winfo_screenheight()
            
       
        x = (screen_width - width) // 2
        y = (screen_height - height) // 2
            
        # Set the window's position
        missing_bul.geometry(f"+{x}+{y}")
        missing_bul.minsize(width, height)
        missing_bul.maxsize(width, screen_height)
        
      
        # Get the directory of the script
        if getattr(sys, 'frozen', False):
        # Running as a bundled executable
            script_directory = sys._MEIPASS
        else:
            # Running as a script
            script_directory = os.path.abspath(os.path.dirname(__file__))
        
        station_dir = os.path.join(script_directory, 'assets', 'station_names.xlsx')
        df = pd.read_excel(station_dir)
        
        
        icon_path = os.path.join(script_directory,'assets', 'dqr_icon.ico')
        missing_bul.iconbitmap(icon_path)
        
        # Initialize observer to None
        observer = None
            
        # Run the main loop
        missing_bul.mainloop()


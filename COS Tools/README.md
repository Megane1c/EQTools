# Create an executable
To create an executable of this script, use pyinstaller module

### Command
pyinstaller --noconsole --add-data "driver;driver" --add-data "gsheet_credentials;gsheet_credentials" --add-data "icons;icons" --add-data "PhaseData;PhaseData" --add-data "Stations;Stations" --icon="icons/tools_icon.ico" --name Tools Tools.py

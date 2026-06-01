#define MyAppName "CyANTs"
#define MyAppVersion "0.1.0"
#define MyAppPublisher "Laboratory of Systems Pharmacology @ Harvard"

[Setup]
AppId={{A9D99452-36A2-4C65-9F72-3E4B36A6E6E0}
AppName={#MyAppName}
AppVersion={#MyAppVersion}
AppPublisher={#MyAppPublisher}
DefaultDirName={localappdata}\Programs\CyANTs
DefaultGroupName=CyANTs
DisableProgramGroupPage=yes
OutputDir=..\..\..\dist
OutputBaseFilename=CyANTs_Setup
SetupIconFile=..\assets\cyants_icon.ico
Compression=lzma2
SolidCompression=yes
PrivilegesRequired=lowest
ArchitecturesAllowed=x64
ArchitecturesInstallIn64BitMode=x64
WizardStyle=modern
UninstallDisplayName={#MyAppName}
UninstallDisplayIcon={app}\assets\cyants_icon.ico

[Tasks]
Name: "desktopicon"; Description: "Create a desktop shortcut"; GroupDescription: "Additional shortcuts:"; Flags: unchecked

[Files]
Source: "..\..\..\*"; DestDir: "{app}"; Excludes: ".git\*,build\*,dist\*,__pycache__\*,*.pyc,.DS_Store"; Flags: ignoreversion recursesubdirs createallsubdirs
Source: "..\..\..\dist\CyANTs.exe"; DestDir: "{app}"; Flags: ignoreversion skipifsourcedoesntexist
Source: "..\runtime\launch_cyants_gui.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\runtime\cyants_conda.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\runtime\install_or_repair_cyants.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\runtime\update_cyants.bat"; DestDir: "{app}"; Flags: ignoreversion
Source: "..\assets\cyants_icon.ico"; DestDir: "{app}\assets"; Flags: ignoreversion

[Icons]
Name: "{group}\CyANTs GUI"; Filename: "{app}\launch_cyants_gui.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\cyants_icon.ico"
Name: "{group}\Install or Repair CyANTs"; Filename: "{app}\install_or_repair_cyants.bat"; Parameters: """{app}"" cyants"; WorkingDir: "{app}"; IconFilename: "{app}\assets\cyants_icon.ico"
Name: "{group}\Update CyANTs"; Filename: "{app}\update_cyants.bat"; WorkingDir: "{app}"; IconFilename: "{app}\assets\cyants_icon.ico"
Name: "{autodesktop}\CyANTs GUI"; Filename: "{app}\launch_cyants_gui.bat"; WorkingDir: "{app}"; Tasks: desktopicon; IconFilename: "{app}\assets\cyants_icon.ico"

[Run]
Filename: "{app}\install_or_repair_cyants.bat"; Parameters: """{app}"" cyants --no-pause"; Description: "Create or update the CyANTs conda environment now"; Flags: postinstall shellexec waituntilterminated skipifsilent
Filename: "{app}\launch_cyants_gui.bat"; Description: "Launch CyANTs GUI"; Flags: postinstall shellexec nowait skipifsilent

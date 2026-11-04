import os
import subprocess

blacklist = [
    '119732', '106521', '117021', '100610', '117930', '103111', '103010', 
    '115017', '120717', '115724', '118124', '120515', '100307', '105014', 
    '113619', '112516', '120212', '107018', '119025', '111413', '114217', 
    '112819', '116726', '100408', '101107', '118932', '112920', '112112', 
    '118023', '113821', '106319', '115219', '102513', '114318', '115320', 
    '111312', '101309', '109123', '118528', '102109', '101410', '117324', 
    '113922', '119833', '113316', '102715', '108828'
]

local_path = os.path.join(os.getcwd(), "original_tractograms")
remote_path = os.environ["REMOTE_TRACTOGRAMS_PATH"]  # e.g. user@host:/path/to/original_tractograms/
password = os.environ["REMOTE_SSH_PASSWORD"]  # never hard-code credentials

subs = os.listdir(local_path)

for sub in subs:
    if sub not in blacklist:
        local_sub_path = os.path.join(local_path, sub)
        command = [
            "sshpass", "-p", password, "scp", "-r", local_sub_path, remote_path
        ]
        
        # Run the command
        print("se")
        result = subprocess.run(command, capture_output=True)
        print("relargoooooo")
        # Print result
        if result.returncode == 0:
            print(f"Successfully copied {sub} to the remote server.")
        else:
            print(f"Failed to copy {sub}. Error: {result.stderr.decode('utf-8')}")

import subprocess

try:
    subprocess.check_output('nvidia-smi')
    print("yes NVIDIA GPU")
except Exception:
    print("no NVIDIA GPU")
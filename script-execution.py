import os
import subprocess


from pathlib import Path


folders = ["notified-alto-santo", "notified-limoeiro"]
algorithms = ['e']
commands = []
folders = [
            "results-exp-default",
            "results-exp-frac-cut",
            "results-exp-default-warm-s",
            "results-exp-frac-cut-warm-s",
            "results-mtz",
        ]

for direc in folders:
    Path(direc).mkdir(parents=True, exist_ok=True)
    
for f in folders:
    instance = os.listdir(f)
    for inst in instance:
        for alg in algorithms:
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-default/" + inst + " " + alg + " 120 8000 1 0 0 0")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-frac-cut/" + inst + " " + alg + " 120 8000 1 0 0 1")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-default-warm-s/" + inst + " " + alg + " 120 8000 1 0 1 0")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-frac-cut-warm-s/" + inst + " " + alg + " 120 8000 1 0 1 1")
            commands.append("./dparp " + f + "/" + inst + " " + "results-mtz/" + inst + " c 120 8000 1 0 0 0")

folders = ["campinas-random"]       
for f in folders:
    instance = os.listdir(f)
    for inst in instance:
        for alg in algorithms:
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-default/" + inst + " " + alg + " 120 8000 1 1 0 0")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-frac-cut/" + inst + " " + alg + " 120 8000 1 1 0 1")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-default-warm-s/" + inst + " " + alg + " 120 8000 1 1 1 0")
            commands.append("./dparp " + f + "/" + inst + " " + "results-exp-frac-cut-warm-s/" + inst + " " + alg + " 120 8000 1 1 1 1")
            commands.append("./dparp " + f + "/" + inst + " " + "results-mtz/" + inst + " c 120 8000 1 1 0 0")
            
for c in commands:
    print(c)
    p = subprocess.Popen(c, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    msg, err = p.communicate()
    if msg:
        print(msg)
    print("OK!!")

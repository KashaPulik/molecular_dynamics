import matplotlib.pyplot as plt

input_file = "./data/data8000.txt"
output_file = "pgraph.jpg"

processes = []
times = []

with open(input_file, 'r') as file:
    for i, line in enumerate(file):
        time, proc = map(float, line.split())
        if i < 5:
            times.append(time)
            processes.append(int(proc)) 

initial_time = times[0]
odin = [1, 4, 8, 12, 16, 20, 24, 28, 32]
more = range(0, 36, 4)[1:]
speedup = [initial_time / t for t in times]

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

plt.figure(figsize=(10, 10))
plt.tick_params(axis='both', width=4)

labels = ['(1x1)', '(4x1)', '(4x2)', '(4x4)', '(4x8)']
plt.xticks(processes, labels)
plt.yticks(odin)
plt.xticks(odin)

plt.plot(processes, speedup, marker='o', linestyle='-', color='r', label="N₂ = 8000")

input_file = "./data/data3375.txt"

processes = []
times = []

with open(input_file, 'r') as file:
    for i, line in enumerate(file):
        time, proc = map(float, line.split())
        if i < 5:
            times.append(time)
            processes.append(int(proc)) 

initial_time = times[0]
speedup = [initial_time / t for t in times]

plt.plot(processes, speedup, marker='v', linestyle='-', color='g', label="N₁ = 3375")

plt.plot(range(1,33),range(1,33),'-',c="blue",linewidth=0.5,label="Линейное ускорение")


plt.xlabel('(Количество узлов)x(Количество процессов)', fontsize=14)
plt.ylabel('Ускорение', fontsize=14)

plt.grid(True, fillstyle='full')
plt.legend()

plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
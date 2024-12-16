import matplotlib.pyplot as plt

input_file = "data.txt"
output_file = "graph1.jpg"

n = []
times = []

with open(input_file, 'r') as file:
    for i, line in enumerate(file):
        time, proc = map(float, line.split())
        if i < 11:
            times.append(time)
            n.append(int(proc)) 

plt.figure(figsize=(10, 10))

plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']

plt.xlabel('Количество частиц n', fontsize=14)
plt.ylabel('Время, с', fontsize=14)

plt.plot(n, times, marker='o', linestyle='-', color='b', label='O(n²)')

plt.grid(True, fillstyle='full')
plt.legend()

plt.savefig(output_file)

print(f"График сохранён в файл {output_file}")
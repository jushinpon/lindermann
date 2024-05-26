import pandas as pd

# 加載數據
file_path = '/mnt/data/Temp_linder_lindemann.csv'
data = pd.read_csv(file_path)

# 重命名列名
data.columns = ['Temperature', 'Lindemann Parameter']

# 設定Lindemann準則臨界值
critical_value = 0.1

# 找到首次達到臨界值的溫度
melting_point = data[data['Lindemann Parameter'] >= critical_value].iloc[0]['Temperature']

print(f"Estimated Melting Point: {melting_point} K")

# 繪圖
import matplotlib.pyplot as plt

plt.figure(figsize=(10, 6))
plt.plot(data['Temperature'], data['Lindemann Parameter'], marker='o', linestyle='-', color='b')
plt.axhline(y=critical_value, color='r', linestyle='--', label=f'Critical Value ({critical_value})')
plt.axvline(x=melting_point, color='g', linestyle='--', label=f'Melting Point ({melting_point} K)')
plt.xlabel('Temperature (K)')
plt.ylabel('Lindemann Parameter')
plt.title('Temperature vs Lindemann Parameter')
plt.legend()
plt.grid(True)
plt.show()

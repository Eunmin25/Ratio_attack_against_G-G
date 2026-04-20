import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# CSV 파일 불러오기
df = pd.read_csv('correlation_data.csv')

# mean_ratio 히스토그램
plt.figure(figsize=(8,4))
sns.histplot(df['mean_ratio'], bins=50, kde=True)
plt.title('Mean Ratio Histogram')
plt.xlabel('mean_ratio')
plt.ylabel('Frequency')
plt.tight_layout()
plt.savefig('mean_ratio_histogram.png')
plt.close()

# mean_ratio boxplot (s_ij별로)
plt.figure(figsize=(8,4))
sns.boxplot(x='s_ij', y='mean_ratio', data=df)
plt.title('Mean Ratio Boxplot by s_ij')
plt.tight_layout()
plt.savefig('mean_ratio_boxplot.png')
plt.close()

print('시각화 결과: mean_ratio_histogram.png, mean_ratio_boxplot.png 파일로 저장됨.')

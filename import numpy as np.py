import numpy as np

# Sigma_full.log에서 행렬만 추출 (CSV 형식)
sigma = np.loadtxt('Sigma_full.log', delimiter=',', skiprows=1)
# 고유값 계산
eigvals = np.linalg.eigvalsh(sigma)
print(eigvals)
print("최소 고유값:", np.min(eigvals))
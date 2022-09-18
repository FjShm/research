# 系の生成手順

1. 初期配置をランダムに作成
2. セルサイズ調整 && エネルギー最小化
3. EQ, NVT, NPT, PROD(NVT)を繰り返してbonded potentials計算
4. セルサイズ調整(密度がtargetになるように)
5. EQ, NVT, NVT(long)を繰り返してnon-bonded
6. NVTでpressure correction

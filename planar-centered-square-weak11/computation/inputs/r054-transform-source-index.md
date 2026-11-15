# R54 左反射最小输入来源索引

- 原始输入：`computation/inputs/r041_final_0036.json`；44767 bytes；SHA-256 `ab714c7f457ccf69f68fd629c947021faab16631c8995caa16b67c1ca2752db3`。原 R41 输入的项目内既有副本，不重制。
- 反射输入（待主代理整合）：`computation/inputs/r054_r041_left_reflected.json`；44845 bytes；SHA-256 `deb65c6011d01f61253cb1db6b4b0e919f014ad32711677b9ee3fd80645694bf`；来源工作区 `pcw-astra-reflection`。
- 原变换记录：`computation/results/double-reflection/r041-transform.json`；1097 bytes；SHA-256 `beb0fcd72ef3ed181ef14ab718cc3c828738ca8f129e35c14a973c3590d9d766`；来源工作区 `pcw-astra-reflection`。其中绝对路径记录原执行环境，不应当作为整合后路径。
- 本次独立只读核验脚本：`computation/verify_r054_transform.py`；SHA-256 `8ef39c2e29137ba12c9434b27447b125d436a687f32fe37878e08658a9b8850b`。Python 3.15.0，退出0；204条按原顺序逐点满足 x'=X−x、y 与质量不变，λ、Px、Py、总质量不变；X=979821360288625294630001/150000000000000000000000。仅验证输入变换，不验证接缝、halo 或常数。

整合后从 topic 根目录复现（只读取两个输入，打印到标准输出）：

```sh
python3 -B computation/verify_r054_transform.py computation/inputs/r041_final_0036.json computation/inputs/r054_r041_left_reflected.json
```

原水平列索引 i≥0 在反射后为 i'≤0，恒等式 X−(α+iPx)=(X−α)+(−i)Px；未做 modulo。冻结 R46 材料沿用既有副本，不在本索引中复制。

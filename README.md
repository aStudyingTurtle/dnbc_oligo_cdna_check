# oligo_cdna_check

DNBelab 单样本 oligo—cDNA barcode 配套审查工具。输入各一个 R1 FASTQ（支持 gzip）；多个文件需用户预先合并。

## 构建与运行

需要 C++17、CMake 和 zlib。在项目目录执行：

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 2
./oligo_cdna_check --sample SAMPLE --cdna-r1 cdna.fq.gz \
  --oligo-r1 oligo.fq.gz --reads 1M --out-prefix SAMPLE
```

已提供可执行文件 `./oligo_cdna_check`；重新构建的程序位于 `build/oligo_cdna_check`。

- cDNA 全量构建 top-50,000 barcode 参考；oligo 取前 `1M`（默认）、`10M` 或 `all`，不是随机抽样。
- 匹配率＝命中参考的 oligo reads / 有效 oligo reads；低于 10% 时警告并标记 `LOW_MATCH_RATE`。
- 输出 JSON、TSV；输出父目录须存在，不覆盖已有结果。不修改原始数据。
- 该指标不是 mapping rate，也不能独立证明样本身份。

**GPT6.1 贡献：代码实现、性能优化及文档整理。**

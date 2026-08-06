# 参照論文

**PDF はこのリポジトリに入っていません。** 各自で下記から取得してください
（`bash papers/fetch_papers.sh` で arXiv 版をここに落とせます）。

再配布しない理由:

- 出版社版 (ApJ) は **© AAS / IOP Publishing**。AAS 誌が CC-BY になったのは
  2022 年からで、それ以前の論文は対象外。
- arXiv 版は **arXiv.org perpetual, non-exclusive license 1.0**。これは
  「**arXiv が**配布してよい」という許諾であって、第三者の再配布は含まない。

数値表（Table 2 など）を `docs/01_paper_analysis.md` に転記し、
教材中で比較対象として使っているのは、**事実・データの出典明記つきの引用**です。

---

## 主教材（この講習会が再現する論文）

**Warren, H. P., Winebarger, A. R., & Brooks, D. H. 2012,
"A Systematic Survey of High-Temperature Emission in Solar Active Regions",
ApJ, 759, 141**

- DOI: https://doi.org/10.1088/0004-637X/759/2/141
- arXiv: https://arxiv.org/abs/1204.3220
- NASA ADS: https://ui.adsabs.harvard.edu/abs/2012ApJ...759..141W

Table 1（15 活動領域）と Table 2（region 7 の 22 輝線 + AIA 94 Å の実測値）が
この教材の「答え合わせ」の基準になります。

`scripts/extract_paper_boxes.py` は、この PDF の Figure 1–3 から
inter-moss 箱の座標を実測します。**PDF を手元に置いてから**実行してください:

```bash
python scripts/extract_paper_boxes.py "papers/Warren_2012_....pdf"
```

## 参照論文（`refs/`）

| ファイル名 | 論文 |
|---|---|
| `arxiv_1009.5976.pdf` | Warren, Brooks & Winebarger 2011, ApJ **734**, 90<br>"Constraints on the Heating of High Temperature Active Region Loops: Observations from Hinode and SDO"<br>https://doi.org/10.1088/0004-637X/734/2/90 ・ https://arxiv.org/abs/1009.5976 |
| `arxiv_1106.5057.pdf` | Winebarger, Schmelz, Warren, Saar & Kashyap 2011, ApJ **740**, 2<br>"Using a Differential Emission Measure and Density Measurements in an Active Region Core to Test a Steady Heating Model"<br>https://doi.org/10.1088/0004-637X/740/1/2 ・ https://arxiv.org/abs/1106.5057 |
| `arxiv_1107.4480.pdf` | Tripathi, Klimchuk & Mason 2011, ApJ **740**, 111<br>"Emission Measure Distribution and Heating of Two Active Region Cores"<br>https://doi.org/10.1088/0004-637X/740/2/111 ・ https://arxiv.org/abs/1107.4480 |

教材のどこで使っているか:

- **Warren+2011**: Fe XIII の R = 1.87 / 1.90（我々と論文の食い違いが
  この論文でも同じであることの根拠）、inter-moss の輝線比
- **Winebarger+2011**: inter-moss 箱のサイズ（5″ × 25″ という細長い箱の流儀）
- **Tripathi+2011**: 別の活動領域の inter-moss 領域の輝線比
  （モジュール 5 で「自然なばらつき」仮説を検証するときの比較対象）

## その他の出典

- **CHIANTI** 原子データベース（`work/gofnt_chianti901*.txt` の計算に使用）:
  Dere et al. 1997, A&AS **125**, 149 / Dere et al. 2019, ApJS **241**, 22（v9.0.1）
  https://www.chiantidatabase.org/
- **PINTofALE** (MCMC_DEM): Kashyap & Drake 1998, ApJ **503**, 450
- **demregpy**: Hannah & Kontar 2012, A&A **539**, A146
- **eispac**: Weberg, Warren et al. 2023, JOSS **8**, 4914

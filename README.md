# str_kernel
文字列カーネルのRust実装です．  
[こちらの記事](https://jetbead.hatenablog.com/entry/20130406/1365183300)を参考にさせていただきました．

## カーネル
- LCS（最長共通部分列）
- スペクトルカーネル
- 固定長部分列カーネル
- ギャップ加重部分列カーネル

## Python
### インストール
1. `pip install maturin`. 
2. `maturin build`. 
3. `pip install target/wheels/[file_name].whl`

### Functions
- `gram_matrix_lcs(List[Str], List[Str])`
- `gram_matrix_spectrum(p, List[Str], List[Str])`
- `gram_matrix_psubseq(p, List[Str], List[Str])`
- `gram_matrix_gapsubseq(lambda, p, List[Str], List[Str])`

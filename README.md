# str_kernel
文字列カーネルのRust実装です．  
[こちらの記事](https://jetbead.hatenablog.com/entry/20130406/1365183300)を参考にさせていただきました．

## カーネル
- LCS（最長共通部分列）
- スペクトルカーネル

## Python
### インストール
1. `pip install maturin`. 
2. `maturin build`. 
3. `pip install target/wheels/[file_name].whl`

### Functions
- `gram_matrix_lcs(List[Str], List[Str])`
- `gram_matrix_spectrum(p, List[Str], List[Str])`

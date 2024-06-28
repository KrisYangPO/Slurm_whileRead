#!/bin/bash
# config_readSampleFormula.sh
# version: framework 1

# Annotation
# ==============================================================================
: '
1. 建立 SRR 所有相關資訊，像是 SRR number, Layout, Strategy, genome 等
   把所有資訊都儲存在一個 sample formula 當中，就可以藉由讀取 sample formula，
   判斷這個檔案的後續分析，決定他的參數內容。

2. while read -r var1 var2 可以有多種變化，可以試著用看看。

3. 這將可以用在不只 SRR download，後續 RNA-seq, ChIP-seq 都可以使用，
   甚至是在一開始就能自動一次判斷所有 sample 該做哪些分析，哪些參數跟選項。
'

# While Loop read files:
# 去學習 While Loop !!!
# example:

while IFS=$'\t' read -r col1 col2
do
    # 將第一欄和第二欄分別儲存到變數中
    var1="$col1"
    var2="$col2"

    # 顯示變數值 (可以根據需求替換為其他操作)
    echo "Column 1: $var1"
    echo "Column 2: $var2"
done < file.txt

#

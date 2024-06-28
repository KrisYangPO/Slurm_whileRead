#!/bin/bash
# Version:2


# !!版本修改：
: '
(解決) 1. 用外部變數帶入 PATH 可以設計多個 pipeline 所需要用到的 tools
(解決) 2. 要設計判斷式，如果 tools path 已經在 PATH 裡面要跳過。
3. 嘗試將更新的 PATH 裡的路徑逐行列出。
'

# Annotation 使用說明
# ==============================================================================
: '
**** 特別注意 ****
執行這個檔案要用 source exportPATH.sh
才能將 export PATH 執行結果儲存。

建立好 toolsPATH.txt，
裡面有 #註解tools 的名稱，還有他的絕對路徑 (下一行)。

要將路徑萃取出來，先用 grep 抓取以 / 開頭的行，
用 paths=(grep command) 的方式將 grep 出來的東西變成陣列。
就可以用 paths[@] 遞迴陣列的方式丟給 forloop 迭代，
每次執行迴圈就 export PATH=${PATH}:${p} 新增路徑到 PATH 裡面。

紀錄：
1. 用外部變數帶入 PATH 可以設計多個 pipeline 所需要用到的 tools
2. 利用判斷式 if [[$PATH =~ $p]] 判斷 $PATH 是否包含 $p，如果有就回報。沒有就 export。

'


# 輸入 config_toolsPATH.sh (裡面含有所有 tools 絕對路徑)
toolsPATH=$1

# 顯示舊的 PATH
echo -e "Old PATH content:\n"${PATH}"\n\n"

# 抓取 toolsPATH
paths=($(grep '^/' ${toolsPATH}))


# 迴圈 export tool PATH 到 PATH 當中
# 另外判斷 PATH 是不是已經包含 tools path (p)了。
for p in ${paths[@]}; do

  if [[ $PATH =~ $p ]]; then
    echo -e 'tools PATH:\n'${p}'\nalready present in PATH'
  else
    echo 'PATH has been added Path: '${p}
    export PATH=${PATH}:${p}

  fi
  sleep 0.5s; done


# 顯示更新 PATH
echo -e "\n\nRenewed PATH content: \n"${PATH}"\n\n"

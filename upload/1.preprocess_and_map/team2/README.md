# 统计画图部分
## 工具
使用python进行表格的处理和画图，需要额外安装package: plotly

```
pip install plotly
```
可使用jupyter notebook打开.ipynb文件，逐个代码框运行即可查看表格和绘图。

## 步骤
### 表格处理
首先将summary.txt读入，使用numpy和pandas进行数据的处理，以及计算平均值，重新写入.csv文件。注意平均值的计算要分别统计绝对数值再求比例进而求平均

### 绘图
分别绘制mapping ratio和length distribution的图
#### mapping ratio
这部分尝试了pie chart来展示整体的平均ratio，又尝试使用plotly、matplotlib的animation模块、以及seaborn的boxplot分别展示细节的样本的ratio分布。

#### length distribution
这部分尝试了使用matplotlib的3D模块、plotly以及折线图展示length distribution。

可以发现几种方法各有优劣，对数据进行了不同维度的展示，其中plotly的交互方式非常先进，可以用鼠标查看更多数据的细节。




# Auto Protein Mutation Analyzer (AMPA)
![APMA](Figure/APMA.png)
<!DOCTYPE html>
<html lang="en" style="height: 100%">
<head>
  <meta charset="utf-8">
</head>
<body style="height: 100%; margin: 0">
  <div id="container" style="height: 100%"></div>

  
  <script type="text/javascript" src="https://registry.npmmirror.com/echarts/5.5.0/files/dist/echarts.min.js"></script>

  <script type="text/javascript">
    var dom = document.getElementById('container');
    var myChart = echarts.init(dom, null, {
      renderer: 'canvas',
      useDirtyRect: false
    });
    var app = {};
    
    var option;

    option = {
  graphic: {
    elements: [
      {
        type: 'text',
        right: 10,
        top: 10,
        style: {
          text: 'APMA',
          fontSize: 50,
          fontWeight: 'bold',
          lineDash: [0, 200],
          lineDashOffset: 0,
          fill: 'transparent',
          stroke: '#000',
          lineWidth: 1
        },
        keyframeAnimation: {
          duration: 5000,
          loop: true,
          keyframes: [
            {
              percent: 0.7,
              style: {
                fill: 'transparent',
                lineDashOffset: 200,
                lineDash: [200, 0]
              }
            },
            {
              // Stop for a while.
              percent: 0.8,
              style: {
                fill: 'transparent'
              }
            },
            {
              percent: 1,
              style: {
                fill: 'black'
              }
            }
          ]
        }
      }
    ]
  }
};

    if (option && typeof option === 'object') {
      myChart.setOption(option);
    }

    window.addEventListener('resize', myChart.resize);
  </script>
</body>
</html>


The AMPA is intended to calculate the features of protein mutation including Amino Acid Web features
such as Betweenness Closeness etc. and other features like Effectiveness, Stiffness and Entropy, Co.evolution.

Also, a stacking model is built automatically to detect the best model to classify the categories you provide.

It can also give you basic figures
![Figure_exp](Figure/Figure_exp.png)
The tool can only operate in Windows system.

See [Run.ipynb](./Run.ipynb) to operate the tool

# Install
- Dependency python = 3.10

- To get the tool, run the following code
```
git clone https://github.com/Spencer-JRWang/APMA
```
- To install the dependency, run the following code
```
pip install .
```
- You should install R in you laptop and add R_Home in your environmental parameter
```
R_HOME = your_route_to_R
```
- Install bio3d, igraph and NACEN package in your R
```
install.packages("bio3d")
```
```
install.pacakges("igraph")
```
The NACEN website: http://sysbio.suda.edu.cn/NACEN

# Message

> 📧: jrwangspencer@stu.suda.edu.cn

> Department of Bioinformatics, Medical School of Soochow University

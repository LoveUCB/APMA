<?php
session_start();
$valid_users = array("wangjingran", "*", "root","123456"); // 替换为合法用户的学号
?>

<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <meta name="keywords" content="graph promoting, amino acid Network, Network Proteomics, dynamics network, machine learning">
  <link rel="stylesheet" href="css/style_index.css">
  <link rel="icon" href="figure/web_icon.ico" type="image/x-icon">
  <style>
   @import url("https://fonts.googleapis.com/css2?family=Poppins:wght@200;300;400;500;600;700&display=swap");
   @import url('https://fonts.googleapis.com/css2?family=Noto+Sans:ital,wght@0,400;0,700;1,400;1,700&display=swap');
  * {
        margin: 0;
        padding: 0;
        box-sizing: border-box;
        font-family: "Poppins", sans-serif;
    }
    .menu {
    position: fixed;
    top: 10px;
    right: 10px;
    display: flex;
    gap: 10px;
  }
  .menu a {
    color: #005af0;
    text-decoration: none;
    padding: 10px 15px;
    background-color: #fff;
    border: 1px solid #005af0;
    border-radius: 5px;
    transition: background-color 0.3s, color 0.3s;
  }
  .menu a:hover {
    background-color: #005af0;
    color: #fff;
  }
  .footer {
            position: relative;
            background-color: rgba(138, 38, 31, 0.9);
            color: white;
            padding: 7px;
            text-align: center;
            width: 100%;
            bottom: 0;
        }
.trail { /* className for the trail elements */
    position: absolute;
    height: 6px; width: 6px;
    border-radius: 3px;
    background: teal;
  }
  </style>
</head>


<body style="overflow-y: scroll;">

    <img src="figure/logo-transparent-png.png" width="360" class="imageContainer"/>


  <div id="container_a" style="height: 200%">
    <script type="text/javascript" src="https://registry.npmmirror.com/echarts/5.5.0/files/dist/echarts.min.js"></script>
    <script type="text/javascript">
      var dom = document.getElementById('container_a');
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
          right: 'center',
          top: '50%',
          style: {
            text: 'APMA',
            fontSize: 70,
            fontWeight: 'bold',
            lineDash: [0, 200],
            lineDashOffset: 0,
            fill: 'transparent',
            stroke: '#DA3F34',
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
                  fill: '#DA3F34'
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
  </div>

<h1 style="text-align: center;">AUTO PROTEIN MUTATION ANALYZER</h1>
<div id = trail></div>
<img src="figure/APMA_main.png" width="80%" class="center-image"/>
<br/>
<h1 style="text-align: center;">Prepare for your files</h1>
  <div class="container">
    <div class="left-box">
      <h3>The Protein Structure File (PDB)</h3>
      <p>If you can get your PDB structure with all amino acid sequence from PDB database, use your own structure file.</p>
      <p>Else you can get the structure from Alphafold database.</p>
      <p>The link to Alphafold Database: <a href="https://alphafold.ebi.ac.uk/">https://alphafold.ebi.ac.uk</a> </p>
      <div class="search-container">
          <form onsubmit="return submitForm('alphafold');">
              <label for="searchStringAlphafold">Search Your protein Structure from Alphafold database:</label><br>
              <input type="text" id="searchStringAlphafold" name="searchString" class="search-input" required>
              <button type="submit" class="search-button">Search</button>
          </form>
      </div>
  </div>
    <div class="right-box">
        <h3>TXT format file for mutations</h3>
        <p>Use tabs to separate phenotypes and mutation types. For mutation types, use the amino acid mutation format (e.g. A3I).</p>
        <p>Example: <a href="https://github.com/Spencer-JRWang/APMA/blob/main/data/position.txt">Example data</a></p>
        <p>The link to gnomAD: <a href="https://alphafold.ebi.ac.uk/">https://gnomad.broadinstitute.org</a> </p>
        <div class="search-container">
            <form onsubmit="return submitForm('gnomad');">
                <label for="searchStringGnomad">Search Your gene in gnomAD:</label><br>
                <input type="text" id="searchStringGnomad" name="searchString" class="search-input" required>
                <button type="submit" class="search-button">Search</button>
            </form>
        </div>
    </div>
  </div>

<script>
function submitForm(database) {
    var searchString;
    var url;
    if (database === 'alphafold') {
        searchString = document.getElementById("searchStringAlphafold").value.trim();
        if (searchString !== "") {
            url = "https://www.alphafold.ebi.ac.uk/search/text/" + encodeURIComponent(searchString);
            window.open(url, "_blank");
        }
    } else if (database === 'gnomad') {
        searchString = document.getElementById("searchStringGnomad").value.trim();
        if (searchString !== "") {
            url = "https://gnomad.broadinstitute.org/gene/" + encodeURIComponent(searchString);
            window.open(url, "_blank");
        }
    }
    return false; // Prevent form submission
}





</script>
</div>



<h2 style="text-align: center;">Model Construction</h2>
<dev class="container_b">
  <img src="figure/luxian.png" class="float-image"/>
  <div class="overlay">APMA模型总体构架: 基于现阶段采用生物信息学方法研究共病机制的局限性，我们提出一个整合序列、能量、结构、动力学等多层面的特征参数，从局部和全局的维度综合考虑突变导致不同表型的识别和共病机制的分析方法，采用机器学习算法搭建一个自动的，可解释的，泛化的模型，并将其封装成一个自动化流程，以python包和在线服务器两种形式提供给需要的人使用，为蛋白质共病体系的分子识别和机制探究提供新的研究契机</div>
<h3>
  1)获取同源序列并进行多序列比对
</br>
</br>
  Blastp官网：<a href="https://blast.ncbi.nlm.nih.gov">https://blast.ncbi.nlm.nih.gov</a>
</br>
</br>
  2)对每一个位点进行保守性打分
</br>
</br>
  Rate4site官网：<a href="https://www.tau.ac.il/~itaymay/cp/rate4site.html">https://www.tau.ac.il/~itaymay/cp/rate4site.html</a>
</br>
</br>
  3)进行突变体构建
</br>
</br>  
  FoldX官网：<a href="https://foldxsuite.crg.eu">https://foldxsuite.crg.eu</a>
</br>
</br>
  4)进行氨基酸网络构建
</br>
</br>  
  NACEN官网：<a href="http://sysbio.suda.edu.cn/NACEN">http://sysbio.suda.edu.cn/NACEN</a>
</br>
</br>
  5)计算蛋白质每一个位点的相对可及面积
</br>
</br>
  6)计算蛋白质的弹性网络参数、熵和共演化系数
</br>
</br>  
  Prody官网：<a href="http://prody.csb.pitt.edu">http://prody.csb.pitt.edu</a>
</br>
</br>
  7)机器学习：降维——搜索——打分
</br>
</br>
  8)模型解释：shap可视化
</br>
</br>
SHAP官网：<a href="https://shap.readthedocs.io/en/latest/">https://shap.readthedocs.io/en/latest/</a>
</h3>
</dev>
</br>
<h2 style="text-align: center;" class="clear">Model Output & Explanation</h2>
</br>

<img src="figure/mo.png" width="60%" class="center-image"/>
</br>
</br>
<div class="iframe-container">
  <iframe src="figure/force_plot_ASD_AC.html" width="2000" height="400" frameborder="0"></iframe>
</div>

</body>



<div class="menu">
  <a href="index.php">HOME</a>
  <a href="APMA.php">RUN</a>
  <a href="login.html">LOGIN</a>
  <a href="guide.html">GUIDE</a>
  <a href="entrez.html">ENTREZ</a>
  <?php
  $valid_users = array("wangjingran", "*", "root","123456");
  if (!isset($_SESSION['user_id']) || !in_array($_SESSION['user_id'], $valid_users)) {
    echo '';
  }
  else{
    echo '<a href="webeditor.php">EDIT</a>';
    echo '<a href="logout.php">LOGOUT</a>';
  }
  ?>

  
</div>



<div class="footer">
  <!-- 底部内容 -->
  Department of Bioinformatics,
  Medical School of Soochow University <br>
  Contact: <a href="mailto:spencer-jrwang@foxmail.com">spencer-jrwang@foxmail.com</a><br>
  Source Code: <a href="https://github.com/Spencer-JRWang/APMA">Github</a>
</div>
</html>

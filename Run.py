# -*- coding: utf-8 -*-

"""

@ author: Jingran Wang

@ Email: jrwangspencer@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

@ GitHub: https://github.com/Spencer-JRWang/APMA


"""

#############################################
### Introduction of Run Module
#
# @ This module is to control the main progress
#
#############################################


from APMA import APMA
# from . import __yourdownloadroute__
# APMA_path = __yourdownloadroute__
# sudo chown -R www-data:www-data /home/wangjingran/APMA
import os
import uuid
import shutil
import glob
import traceback
import datetime
import re
import requests
import time
import PyPDF2

# fetch start time
current_datetime = datetime.datetime.now()

# print start time
print("***** APMA Start at: ", current_datetime)


def find_consecutive_numbers(input_string):
    """
    Finds all consecutive numbers in a given string and returns them as a list of strings.

    Parameters:
    input_string (str): The input string in which to find consecutive numbers.

    Returns:
    list: A list of strings representing the consecutive numbers found in the input string.
    """
    pattern = r'\d+'
    matches = re.findall(pattern, input_string)
    return matches


def download_alphafold_structure(uniprot_id, output_path):
    """
    Downloads the protein structure file from the AlphaFold database based on the given Uniprot ID and saves it to the specified output path.

    Parameters:
    uniprot_id (str): The Uniprot ID of the protein for which the structure file is to be downloaded.
    output_path (str): The path (including filename) where the downloaded structure file should be saved.

    Returns:
    None

    Example:
    download_alphafold_structure("P12345", "P12345.pdb")
    """
    print("[INFO] ...Downloading Structure File from AlphaFold database...")
    base_url = "https://alphafold.ebi.ac.uk/files/AF-"
    file_url = f"{base_url}{uniprot_id}-F1-model_v4.pdb"
    response = requests.get(file_url)

    if response.status_code == 200:
        with open(output_path, 'wb') as file:
            file.write(response.content)
        print(f"[INFO] Structure file for {uniprot_id} downloaded successfully.")
    else:
        raise ValueError(f"[INFO] Failed to download structure file for {uniprot_id}. Status code: {response.status_code}")


def delete_files_in_directory(directory):
    """
    Deletes all files and directories within the specified directory.
    
    Parameters:
    directory (str): The path to the directory from which to delete files and directories.
    
    Returns:
    None
    
    This function traverses the specified directory and deletes all files.
    If a subdirectory is encountered, it recursively deletes its contents as well.
    
    Example:
    delete_files_in_directory("/home/user/directory_to_clean")
    """
    # Traverse all files and subdirectories in the specified directory
    for item in os.listdir(directory):
        item_path = os.path.join(directory, item)
        
        # If the item is a file, delete it
        if os.path.isfile(item_path):
            os.remove(item_path)
        
        # If the item is a directory, recursively call this function
        elif os.path.isdir(item_path):
            delete_files_in_directory(item_path)


def create_task_id_directory(base_directory):
    """
    Creates a unique task directory within the specified base directory.

    This function generates a unique UUID for each task, checks if a directory
    with the same name already exists in the base directory, and creates a new
    directory with the UUID as its name if it does not already exist.

    Parameters:
    base_directory (str): The path to the base directory where the task directories will be created.

    Returns:
    str: The path to the newly created task directory.
    None: If an error occurs during the directory creation.
    """
    while True:
        # Generate a unique task ID
        task_id = str(uuid.uuid4())
        print(f"[INFO] Task Successfully created: {task_id}")
        # Create the path for the task directory
        task_directory = os.path.join(base_directory, task_id)
        
        # Check if a directory with the same name already exists
        if not os.path.exists(task_directory):
            try:
                # Create the directory
                os.makedirs(task_directory)
                print(f"[INFO] Task directory created: {task_directory}")
                return task_id, task_directory
            except OSError as e:
                print(f"[ERROR] Error creating task directory: {e}")
                return None
        else:
            # If a directory with the same name exists, generate a new UUID
            continue


def merge_pdfs_from_folder(folder_path, output_pdf):
    """
    Merge all PDF files from a specified folder into a single PDF file.

    Parameters:
    folder_path (str): Path to the folder containing PDF files to merge.
    output_pdf (str): Output file path for the merged PDF.

    Returns:
    None

    Raises:
    ValueError: If no PDF files are found in the specified folder.
    FileNotFoundError: If the specified folder does not exist.
    """

    # Validate folder_path exists
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f"The folder '{folder_path}' does not exist.")

    # List all PDF files in the folder
    pdf_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith('.pdf')]

    # Check if there are any PDF files in the folder
    if not pdf_files:
        raise ValueError(f"No PDF files found in the folder '{folder_path}'.")

    # Create a PdfFileWriter object to hold the combined PDF
    output_pdf_writer = PyPDF2.PdfWriter()

    # Iterate through each PDF file and add its content to output_pdf_writer
    for pdf_file in pdf_files:
        try:
            with open(pdf_file, 'rb') as f:
                input_pdf_reader = PyPDF2.PdfReader(f)
                num_pages = len(input_pdf_reader.pages)
                

                # Add each page of the input PDF to the output PdfFileWriter
                for page_num in range(num_pages):
                    page = input_pdf_reader.pages[page_num]
                    output_pdf_writer.add_page(page)
        except FileNotFoundError:
            raise FileNotFoundError(f"[ERROR] Input PDF file '{pdf_file}' not found.")

    # Write the merged PDF content to a file
    with open(output_pdf, 'wb') as f:
        output_pdf_writer.write(f)

    print(f"[INFO] PDFs merged successfully into {output_pdf}")


def create_inprocess_view_html(task_id, protein_name, email):
    """
    Create inprocess view html file
    This file can tell the status, id, protein name, email etc information
    Use chart.js to preview the user's data

    Parameters:
    - task_id: user's task id
    - protein_name: user's uniprot ID or user's pdb file prefix
    - email: user's email
    """
    # move position.txt to /var/www/html/deePheMut
    shutil.copyfile('/home/wangjingran/APMA/data/position.txt', '/var/www/html/deePheMut/position.txt')
    ## fit the mutation file to chart.js
    file_path = "/var/www/html/deePheMut/position.txt"
    with open(file_path, 'r') as file:
        lines = file.readlines()
    if lines and lines[-1].endswith('\n'):
        lines[-1] = lines[-1].strip()
    with open(file_path, 'w') as file:
        file.writelines(lines)
    
    html_code = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <title>deePheMut progress</title>
    <link rel="icon" href="../../figure/web_icon.ico" type="image/x-icon">
    <link rel="stylesheet" href="../../inprocess.css">
</head>
<body>
<img src="../../figure/logo-transparent-png.png" width="360" class="imageContainerr"/>
<div id="container_a" style="height: 80%">
    <script type="text/javascript" src="https://registry.npmmirror.com/echarts/5.5.0/files/dist/echarts.min.js"></script>
    <script type="text/javascript" src="../../echarts_cartoon.js"></script>
</div>
<h4 style="text-align: center;">deep and precise prediction of mutations to different phenotypes</h4>
<br/>


<div class="container_box">
    <div class="card">
        <div class = "image">
            <canvas id="myChart"></canvas>
            <script src="../../barchart.js"></script>
        </div>
        <div class="content">
            <h1> Status: Running</h1>
            <hr style="height: 3px; background-color: black; border: none;">
            <br>
            <h2> Your Query ID: {task_id}</h2>
            <br>
            <h2> Your Protein is: {protein_name}</h2>
            <br>
            <h2> Your Email is: {email[0]}</h2>
            <br>
        </div>
    </div>
</div>
<div class="menu">
    <a href="../../index.php">HOME</a>
    <a href="../../run.php">RUN</a>
    <a href="../../guide.html">GUIDE</a>
    <a href="../../feedback.php">FEEDBACK</a>
    <button class="btn_change_theme" id="toggleBackgroundBtn">Change Theme</button>
</div>
<script>
    document.getElementById('toggleBackgroundBtn').addEventListener('click', function() {{
        document.body.classList.toggle('dark-theme');
    }});
</script>
<div class="footer">
    Department of Bioinformatics,
    Medical School of Soochow University <br>
    Contact us: <a href="mailto:spencer-jrwang@foxmail.com">spencer-jrwang@foxmail.com</a><br>
    Source Code: <a href="https://github.com/Spencer-JRWang/APMA">Github</a>
</div>
</body>
</html>
"""
    # write to index.php file
    with open(f'/var/www/html/deePheMut/user_data/{task_id}/index.html', 'w') as file:
        file.write(html_code)
    print(f"[INFO] inprocess html file is written")
    return None


def create_error_view_html(task_id, protein_name, email):
    """
    Create inprocess view html file
    This file can tell the status, id, protein name, email etc information
    Use chart.js to preview the user's data

    Parameters:
    - task_id: user's task id
    - protein_name: user's uniprot ID or user's pdb file prefix
    - email: user's email
    """
    
    html_code = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <title>deePheMut progress</title>
    <link rel="icon" href="../../figure/web_icon.ico" type="image/x-icon">
    <link rel="stylesheet" href="../../inprocess.css">
</head>
<body>
<img src="../../figure/logo-transparent-png.png" width="360" class="imageContainerr"/>
<div id="container_a" style="height: 80%">
    <script type="text/javascript" src="https://registry.npmmirror.com/echarts/5.5.0/files/dist/echarts.min.js"></script>
    <script type="text/javascript" src="../../echarts_cartoon.js"></script>
</div>
<h4 style="text-align: center;">deep and precise prediction of mutations to different phenotypes</h4>
<br/>


<div class="container_box">
    <div class="card">
        <div class = "image">
            <canvas id="myChart"></canvas>
            <script src="../../barchart.js"></script>
        </div>
        <div class="content">
            <h1> Status: Error</h1>
            <hr style="height: 3px; background-color: black; border: none;">
            <br>
            <h2> Your Query ID: {task_id}</h2>
            <br>
            <h2> Your Protein is: {protein_name}</h2>
            <br>
            <h2> Your Email is: {email[0]}</h2>
            <br>
        </div>
    </div>
</div>
<div class="menu">
    <a href="../../index.php">HOME</a>
    <a href="../../run.php">RUN</a>
    <a href="../../guide.html">GUIDE</a>
    <a href="../../feedback.php">FEEDBACK</a>
    <button class="btn_change_theme" id="toggleBackgroundBtn">Change Theme</button>
</div>
<script>
    document.getElementById('toggleBackgroundBtn').addEventListener('click', function() {{
        document.body.classList.toggle('dark-theme');
    }});
</script>
<div class="footer">
    Department of Bioinformatics,
    Medical School of Soochow University <br>
    Contact us: <a href="mailto:spencer-jrwang@foxmail.com">spencer-jrwang@foxmail.com</a><br>
    Source Code: <a href="https://github.com/Spencer-JRWang/APMA">Github</a>
</div>
</body>
</html>
"""
    # write to index.php file
    with open(f'/var/www/html/deePheMut/user_data/{task_id}/index.html', 'w') as file:
        file.write(html_code)
    print(f"[INFO] inprocess html file is written")
    return None






def final_view_index_html(task_id):
    """
    Create output preview html file
    This file can tell the status, id, protein name, email etc information
    Use chart.js to preview the user's data

    Parameters:
    - task_id: user's task id
    - protein_name: user's uniprot ID or user's pdb file prefix
    - email: user's email
    """

    # Copy Outcome Folder to task folder
    outcome_folder = '/home/wangjingran/APMA/Outcome'
    destination_folder = f'/var/www/html/deePheMut/user_data/{task_id}'
    folder_name = 'Outcome'
    destination_path = os.path.join(destination_folder, folder_name)
    shutil.copytree(outcome_folder, destination_path)

    # Copy Outcome.zip to the task folder
    shutil.copy('/home/wangjingran/APMA/Email/Outcome.zip', destination_folder)

    def remove_suffix(text, suffix):
        if text.endswith(suffix):
            return text[:-len(suffix)]
        return text

    def remove_prefix(text, prefix):
        if text.startswith(prefix):
            return text[len(prefix):]
        return text

    # define iframe code
    iframe_code = ""
    search_pattern = os.path.join(f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/Explain', '*.html')
    html_files = glob.glob(search_pattern)
    print(html_files)
    html_file_names = [os.path.basename(file) for file in html_files]
    print(html_file_names)
    for html_file in html_file_names:
        iframe_code += f"""<p>{remove_prefix(remove_suffix(html_file,'.html'), 'LightGBM_force_plot_')}</p>\n
        <iframe src='Outcome/Figure/Explain/{html_file}' title='iframe_force'></iframe>\n"""
    
    # define score boxplot code
    box_plot_code = ""
    search_pattern = os.path.join(f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Score', '*.txt')
    score_files = glob.glob(search_pattern)
    score_names = [os.path.basename(file) for file in score_files]
    for score_name in score_names:
        box_plot_code += f"{{ file: './Outcome/Score/{score_name}', title: '{remove_prefix(remove_suffix(score_name,'.txt'), 'LightGBM_')}' }},\n"
    html_code = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <link rel="icon" href="../../figure/web_icon.ico" type="image/x-icon">
    <script src="../../plotly-latest.js"></script>
    <title>Task Preview</title>
    <link rel="stylesheet" href="../../outcome_preview.css">
</head>
<body>
<h1>Task successfully finished</h1>

<div id="pdf-container">
    <h2>Feature Calculation</h2>
    <div class="pdf-item">
        <embed src="./Outcome/Figure/combined_plots.pdf" type="application/pdf">
        <p><strong>Figure 1| Mutation parameter distribution plot and PCA plot</strong>,
            the first 15 images represent the performance of phenotypes in different parameters.
            When the sample size is greater than 30, t test in statistical test will be adopted,
            otherwise wilcoxon test will be used, where ** represents 0.01 < pvalue < 0.05, *** indicates pvalue < 0.001,
            the last picture shows the results of principal component analysis of mutations of different phenotypes.
            The distribution of principal component analysis shows the distribution of different phenotypic mutations on the first two principal components.</p>
        <button onclick="downloadFile('./Outcome/Figure/combined_plots.pdf')">Download Figure 1</button>
        <button onclick="downloadFile('./Outcome/paras.txt')">Download parameters txt</button>
    </div>

    <div class="pdf-item">
        <embed src="./Outcome/Figure/dynamic.pdf" type="application/pdf">
        <p><strong>Figure 2| Destribution of dynamic network based on ANM/ENM model</strong>,
            the plot </p>
        <button onclick="downloadFile('./Outcome/Figure/dynamic.pdf')">Download Figure 1</button>
    </div>

    <div class="pdf-item">
        <embed src="./Outcome/Figure/ROC/Feature/feature_roc.pdf" type="application/pdf">
        <p><strong>Figure 2| Mutation parameter distribution plot and PCA plot</strong>,
            the first 15 images represent the performance of phenotypes in different parameters.
            When the sample size is greater than 30, t test in statistical test will be adopted,
            otherwise wilcoxon test will be used, where ** represents 0.01 < pvalue < 0.05, *** indicates pvalue < 0.001,
            the last picture shows the results of principal component analysis of mutations of different phenotypes.
            The distribution of principal component analysis shows the distribution of different phenotypic mutations on the first two principal components.</p>
        <button onclick="downloadFile('./Outcome/Figure/ROC/Feature/feature_roc.pdf')">Download Figure 2</button>
    </div>

    <h2>Model Construction</h2>

    <div class="pdf-item">
        <embed src="./Outcome/Figure/rfe.pdf" type="application/pdf">
        <p><strong>Figure 3| </strong></p>
        <button onclick="downloadFile('./Outcome/Figure/rfe.pdf')">Download Figure 2</button>
    </div>

    <div class="pdf-item">
        <embed src="Outcome/Figure/spearman.pdf" type="application/pdf">
        <p><strong>Figure 4| Spearman correlation heatmap</strong>, this heatmap illustrates Spearman correlation coefficients computed between protein mutations
            and 15 selected features. Color intensity reflects the strength of correlation: bluer shades indicate positive correlations,
            while redder shades indicate negative correlations. Numeric values adjacent to the color bar denote the magnitude of
            correlation coefficients.</p>
        <button onclick="downloadFile('./Outcome/Figure/spearman.pdf')">Download Figure 2</button>
    </div>

    <div class="pdf-item">
        <pre id="txtContent"></pre>
        <p><strong>Text 1| Model construction log</strong>, this heatmap illustrates Spearman correlation coefficients computed between protein mutations
            and 15 selected features. Color intensity reflects the strength of correlation: bluer shades indicate positive correlations,
            while redder shades indicate negative correlations. Numeric values adjacent to the color bar denote the magnitude of
            correlation coefficients.</p>
        <button onclick="downloadFile('./Outcome/Feature_selection.txt')">Download TXT 1</button>
    </div>

    <h2>Model Explanation</h2>
    <div class="pdf-item">
        <embed src="Outcome/Figure/Explain/shap_summary_plot.pdf" type="application/pdf">
        <p><strong>Figure 5| Machine learning model explanation force plot</strong>, This Force Plot generated using SHAP visualizes the impact of
            features on the model's predictions.  Each feature's contribution to the prediction is represented by the length and direction(colors).
            Positive and negative contributions are indicated by pointing upwards(red) and downwards(blue), respectively. </p>
        <button onclick="downloadFile('./Outcome/Figure/Explain/shap_summary_plot.pdf')">Download Figure 2</button>
    </div>

    <div class="pdf-item">
        {iframe_code}
        <p><strong>Figure 6| Machine learning model explanation force plot</strong>, This Force Plot generated using SHAP visualizes the impact of
            features on the model's predictions.  Each feature's contribution to the prediction is represented by the length and direction(colors).
            Positive and negative contributions are indicated by pointing upwards(red) and downwards(blue), respectively. </p>
    </div>
    <h2>PhenoScores</h2>
    <div class="pdf-item">
        <div id="boxplots-container" class="boxplot-container"></div>
        <script>
            async function loadDataAndDraw(containerId, dataFile, title) {{
                try {{
                    const response = await fetch(dataFile);
                    const data = await response.text();

                    const rows = data.trim().split('\\n').map(row => row.split('\\t'));
                    drawBoxPlots(containerId, rows, title); // 传入标题参数
                }} catch (error) {{
                    console.error('Error fetching or parsing data:', error);
                }}
            }}

            function drawBoxPlots(containerId, data, title) {{
                var categories = {{}};
                var categoryData = [];
                var colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
                    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'];

                data.forEach((row, index) => {{
                    if (index === 0) return; // Skip header row
                    var category = row[0]; // Assuming first column is category
                    var score = parseFloat(row[2]); // Assuming third column is score
                    var mutation = row[3]; // Assuming fourth column is mutation name

                    if (!categories[category]) {{
                        categories[category] = {{
                            scores: [],
                            mutations: [],
                            color: colors[Object.keys(categories).length % colors.length]
                        }};
                    }}
                    categories[category].scores.push(score);
                    categories[category].mutations.push(mutation);
                }});

                for (var category in categories) {{
                    if (categories.hasOwnProperty(category)) {{
                        categoryData.push({{
                            y: categories[category].scores,
                            type: 'box',
                            name: 'Category ' + category,
                            boxpoints: 'all',
                            jitter: 0.3,
                            pointpos: -1.8,
                            marker: {{
                                color: categories[category].color
                            }},
                            line: {{
                                color: categories[category].color
                            }},
                            text: categories[category].mutations,
                            hoverinfo: 'text+y'
                        }});
                    }}
                }}

                var layout = {{
                    title: title, // 使用传入的标题
                    yaxis: {{
                        title: 'Scores'
                    }}
                }};

                Plotly.newPlot(containerId, categoryData, layout);
            }}

            window.addEventListener('DOMContentLoaded', async () => {{
                const dataFiles = [
                    {box_plot_code}
                ];

                const container = document.getElementById('boxplots-container');
                dataFiles.forEach((dataInfo, index) => {{
                    const containerId = `boxplots-${{index + 1}}`;
                    const newDiv = document.createElement('div');
                    newDiv.id = containerId;
                    newDiv.classList.add('boxplot-container');
                    container.appendChild(newDiv);

                    loadDataAndDraw(containerId, dataInfo.file, dataInfo.title); // 传入标题参数
                }});
            }});
        </script>

        <p><strong>Figure 7| </strong></p>
    </div>

    <div class="pdf-item">
        <embed src="./Outcome/Figure/ROC/ML/ml_roc.pdf" type="application/pdf">
        <p><strong>Figure 8| Machine learning model explanation force plot</strong>, This Force Plot generated using SHAP visualizes the impact of
            features on the model's predictions.  Each feature's contribution to the prediction is represented by the length and direction(colors).
            Positive and negative contributions are indicated by pointing upwards(red) and downwards(blue), respectively. </p>
        <button onclick="downloadFile('./Outcome/Figure/ROC/ML/ml_roc.pdf')">Download Figure 2</button>
    </div>
</div>
<button id="download-all" onclick="downloadAllFiles()">Download All Outcomes</button>

<!-- 右侧导航栏 -->
<div id="sidebar">
    <strong>Navigation</strong>
    <ul id="nav-links">
        <!-- JavaScript 会动态生成导航链接 -->
    </ul>
</div>

<script>
    function downloadFile(fileSrc) {{
        const link = document.createElement('a');
        link.href = fileSrc;
        link.download = fileSrc.split('/').pop();
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }}

    function downloadAllFiles() {{
        const link = document.createElement('a');
        link.href = './Outcome.zip'; // 预先存在的ZIP文件路径
        link.download = './Outcome.zip';
        document.body.appendChild(link);
        link.click();
        document.body.removeChild(link);
    }}

    window.onload = function() {{
        // 动态生成导航链接
        const navLinks = document.getElementById('nav-links');
        const headings = document.querySelectorAll('#pdf-container h2');
        headings.forEach((heading, index) => {{
            const id = 'heading-' + index;
            heading.id = id;
            const li = document.createElement('li');
            const a = document.createElement('a');
            a.href = '#' + id;
            a.innerText = heading.innerText;
            li.appendChild(a);
            navLinks.appendChild(li);
        }});

        // 加载文本文件内容
        fetch('./Outcome/Feature_selection.txt')
            .then(response => response.text())
            .then(data => {{
                document.getElementById('txtContent').textContent = data;
            }})
            .catch(error => {{
                console.error('Error fetching the file:', error);
            }});
    }};
</script>
</body>
</html>
"""
    
    # write to index.html file
    with open(f'/var/www/html/deePheMut/user_data/{task_id}/index.html', 'w') as file:
        file.write(html_code)
    print(f"[INFO] outcome HTML file is written")
    return None




# Define Amino Acid list
Amino_acids_list = [
    'A',  # Alanine
    'R',  # Arginine
    'N',  # Asparagine
    'D',  # Aspartic acid
    'C',  # Cysteine
    'Q',  # Glutamine
    'E',  # Glutamic acid
    'G',  # Glycine
    'H',  # Histidine
    'I',  # Isoleucine
    'L',  # Leucine
    'K',  # Lysine
    'M',  # Methionine
    'F',  # Phenylalanine
    'P',  # Proline
    'S',  # Serine
    'T',  # Threonine
    'W',  # Tryptophan
    'Y',  # Tyrosine
    'V',  # Valine
    'X'   # Unknown
]

# The main console
if __name__ == "__main__":
    try:
        # fetch email
        email_list = []
        f = open("/home/wangjingran/APMA/data/email.txt")
        lines = f.readlines()
        for i in lines:
            line = i.strip("\n")
            email_list.append(line)
        f.close()

        # check whether need to download structure from AlphaFold database
        uni_file_path = "/home/wangjingran/APMA/data/uniprot_ID.txt"
        if os.path.exists(uni_file_path):
            submit = False
        else:
            submit = True

        # Get the Uniprot ID from the user
        if submit:
            pass
        else:
            f = open('/home/wangjingran/APMA/data/uniprot_ID.txt', 'r')
            uniprot_id = f.readlines()
            count_id = 0
            for i in uniprot_id:
                if i == "\n":
                    pass
                else:
                    count_id += 1
            if count_id == 1:
                uniprot_id = uniprot_id[0].strip('\n')
            else:
                raise ValueError("[ERROR] We do not support multiple uniprot IDs")
            f.close()
            # Download the structure file
            download_alphafold_structure(uniprot_id, f"/home/wangjingran/APMA/data/{uniprot_id}.pdb")

        # pocess position file
        dict_res_type = {}
        f = open("/home/wangjingran/APMA/data/position.txt","r")
        all = f.readlines()
        all_new = []
        # print(all)
        for i in all:
            if i == '\n':
                pass
            else:
                string_split = i.split("\t")
                string_a = string_split[0]
                string_b = string_split[1]

                # check amino acid
                ori_amino_acid = string_b[0]
                mut_amino_acid = string_b[-2]
                if (ori_amino_acid not in Amino_acids_list) or (mut_amino_acid not in Amino_acids_list):
                    print("[INFO] Original Amino Acid: " + str(ori_amino_acid))
                    print("[INFO] Mutated Amino Acid: " + str(mut_amino_acid))
                    raise ValueError("[ERROR] Unsupported amino acid abbreviation")
                
                # change to FoldX forat
                string_b = string_b[:1] + "A" + string_b[1:]

                # add into dict
                dict_res_type[int(find_consecutive_numbers(string_b)[0])] = string_b[0]
                FoldX_type_string = string_a + "\t" + string_b
                all_new.append(FoldX_type_string)
        f.close()


        f = open("/home/wangjingran/APMA/data/position.txt","w")
        for i in all_new:
            f.write(i)
        f.close()

        del all
        del all_new

        # fetch pdb file
        def print_pdb_files(folder_path):
            # 使用 glob 模块列出文件夹中所有的 .pdb 文件
            pdb_files = glob.glob(os.path.join(folder_path, '*.pdb'))
            for i in pdb_files:
                user_pdb_file = i
            return user_pdb_file

        folder_path = '/home/wangjingran/APMA/data'
        user_pdb = print_pdb_files(folder_path)
        user_protein_name = user_pdb.split("/")[-1].rstrip(".pdb")
        # print(dict_res_type)

        # check position file
        from Feature_Cal.Blast_MSA import extract_sequence_from_pdb
        pdb_seq = extract_sequence_from_pdb(user_pdb)
        pdb_seq = ''.join(pdb_seq)
        for (position, residue) in dict_res_type.items():
            if pdb_seq[position - 1] != residue:
                raise ValueError(f"[ERROR] Position {position} in the position file does not match the PDB sequence.")
        print("[INFO] Pass checking")
        
        # check email
        for i in email_list:
            if '@' not in i:
                raise ValueError("[ERROR] Invalid email address.")

        # After passing several checks
        # Create the task ID and task directory
        # Record User's task and send item ID
        task_id, task_folder = create_task_id_directory('/var/www/html/deePheMut/user_data')

        # create inprocess index.html
        create_inprocess_view_html(task_id, user_protein_name, email_list)

        # Send start email
        from Email.send import send_start_email
        send_start_email(task_id, email_list)

        #####################################################
        # Start the main process
        APMA(
        Protein_name = user_protein_name,
        file_path = "/home/wangjingran/APMA/data/position.txt",
        WT_PDB = user_pdb
        )
        #####################################################

        #####################################################
        # Plot basic figures
        from ML.figure import plot_roc_for_disease_pairs
        plot_roc_for_disease_pairs("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure/ROC/Feature")

        from ML.figure import plot_box
        plot_box("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure")

        from ML.figure import plot_spearman
        plot_spearman("/home/wangjingran/APMA/data/paras.txt","/home/wangjingran/APMA/Outcome/Figure")

        from ML.figure import plot_dynamic_network
        plot_dynamic_network('/home/wangjingran/APMA/data/all_dyn_data.txt', '/home/wangjingran/APMA/data/paras.txt', '/home/wangjingran/APMA/Outcome/Figure/dynamic.pdf')
        #####################################################

        from Email.zip import zip_folder
        zip_folder('/home/wangjingran/APMA/Outcome','/home/wangjingran/APMA/Email/Outcome.zip')
        
        # construct index.html file
        final_view_index_html(task_id)

        # merge pdf files
        merge_pdfs_from_folder(f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/Explain', 
                               f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/Explain/shap_summary_plot.pdf')
        merge_pdfs_from_folder(f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/ROC/Feature', 
                               f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/ROC/Feature/feature_roc.pdf')
        merge_pdfs_from_folder(f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/ROC/ML', 
                               f'/var/www/html/deePheMut/user_data/{task_id}/Outcome/Figure/ROC/ML/ml_roc.pdf')
        
        # print end time
        print("[INFO] APMA ends at: ", current_datetime)

        from Email.send import send_email
        send_email(task_id, email_list)
        

    except Exception as e:
        print(str(e))

        traceback_info = traceback.format_exc()
        print(traceback_info)

        current_datetime = datetime.datetime.now()
        print("[INFO] APMA ends at: ", current_datetime)

        time.sleep(3)

        # create inprocess index.html
        try:
            create_error_view_html(task_id, user_protein_name, email_list)
        except Exception as ee:
            print("[INFO] " + str(ee))
        
        from Email.send import send_error_email
        send_error_email(email_list)


    delete_files_in_directory("/home/wangjingran/APMA/Outcome")
    delete_files_in_directory("/home/wangjingran/APMA/data")

    folder_path = '/home/wangjingran/APMA/FoldX'
    files = os.listdir(folder_path)

    for file_name in files:
        
        if file_name != 'foldx4' and file_name != 'rotabase.txt' and file_name != 'foldx5' and file_name != 'molecules':
            file_path = os.path.join(folder_path, file_name)

            os.remove(file_path)


<?php
session_start();
$valid_users = array("wangjingran", "*", "root","123456"); // 替换为合法用户的学号

if (!isset($_SESSION['user_id']) || !in_array($_SESSION['user_id'], $valid_users)) {
    header('Location: login.html');
    exit();
}
?>


<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>APMA Run Module</title>
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
        html, body {
            height: 100%;
            margin: 0;
        }
        body {
            background-image: url('figure/5uvJN.png');
            background-size: cover;
            background-repeat: no-repeat;
            background-attachment: fixed;
            transition: background-image 0.5s ease-in-out;
        }
        .container {
            height: 110vh;
            width: 100%;
            align-items: center;
            display: flex;
            justify-content: center;
            background-color: transparent;
        }
        .card {
            border-radius: 10px;
            box-shadow: 0 5px 10px 0 rgba(0, 0, 0, 0.3);
            width: 800px;
            height: 1000px;
            background-color: #ffffff;
            padding: 10px 30px 40px;
            transition: background-color 0.5s ease-in-out;
        }
        .card h3 {
            font-size: 22px;
            font-weight: 600;
        }
        .drop_box {
            margin: 10px 0;
            padding: 30px;
            display: flex;
            align-items: center;
            justify-content: center;
            flex-direction: column;
            border: 3px dotted #a3a3a3;
            border-radius: 5px;
        }
        .drop_box_submit {
            display: flex;
            align-items: center;
            justify-content: center;
            flex-direction: column;
            border: none;
            border-radius: 5px;
        }
        .drop_box h4 {
            font-size: 16px;
            font-weight: 400;
            color: #2e2e2e;
        }
        .drop_box p {
            margin-top: 10px;
            margin-bottom: 20px;
            font-size: 12px;
            color: #a3a3a3;
        }
        .btn {
            text-decoration: none;
            background-color: #005af0;
            color: #ffffff;
            padding: 10px 20px;
            border: none;
            outline: none;
            transition: 0.3s;
        }
        .btn:hover {
            text-decoration: none;
            background-color: #ffffff;
            color: #005af0;
            padding: 10px 20px;
            border: none;
            outline: 1px solid #010101;
        }
        .btn_change_theme {
            position: absolute;
            top: 50px;
            right: 10px;
            text-decoration: none;
            background-color: #005af0;
            color: #ffffff;
            padding: 10px 20px;
            border: none;
            outline: none;
            transition: 0.3s;
        }
        .btn_change_theme:hover {
            text-decoration: none;
            background-color: #ffffff;
            color: #005af0;
            padding: 10px 20px;
            border: none;
            outline: 1px solid #010101;
        }
        .btn_submit {
            text-decoration: none;
            background-color: rgba(138, 38, 31, 0.9);
            color: #ffffff;
            padding: 10px 50px;
            border: none;
            outline: none;
            transition: 0.3s;
            font-size: 12px;
        }
        .btn_submit:hover {
            text-decoration: none;
            background-color: #ffffff;
            color: rgba(138, 38, 31, 0.9);
            padding: 10px 20px;
            border: none;
            outline: 1px solid #010101;
        }
        .form input {
            margin: 10px 0;
            width: 100%;
            background-color: #e2e2e2;
            border: none;
            outline: none;
            padding: 12px 20px;
            border-radius: 4px;
        }
        .search-input {
            width: 100%;
            padding: 8px;
            border: 1px solid #ccc;
            border-radius: 4px;
            outline: none;
        }
        .mutation-input {
            width: 100%;
            padding: 20px;
            border: 1px solid #ccc;
            border-radius: 4px;
            outline: none;
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
        .imageContainer {
            position: relative;
            transition: transform 3s ease;
            top: 0;
            left: 0;
            background-color: transparent;
        }
        .imageContainer:hover {
            transform: rotate(360deg);
        }
        a {
            text-decoration: none;
            color: #5dabff;
            transition: color 0.3s;
        }
        a:hover {
            color: #0056b3;
        }
        .clock {
            font-size: 15px;
            font-weight: 600;
            margin-bottom: 15px;
        }

        .dark-theme {
            background-image: url('figure/6193479.jpg');
        }
        .dark-theme .card {
            background-color: #444;
            color: #fff;
        }
        .dark-theme .btn {
            background-color: #555;
            color: #fff;
        }
        .dark-theme .btn:hover {
            background-color: #666;
            color: #fff;
        }
        .dark-theme .btn_change_theme {
            background-color: #555;
            color: #fff;
        }
        .dark-theme .btn_change_theme:hover {
            background-color: #666;
            color: #fff;
        }
        .dark-theme .footer {
            background-color: #222;
        }
        .dark-theme .drop_box h4 {
            color: #fff;
        }
        .dark-theme .drop_box h3 {
            color: #fff;
        }
        .menu {
            position: absolute;
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
    </style>
</head>
<body onload="fetchStatus()">

<div class="menu">
    <a href="index.php">HOME</a>
    <a href="APMA.php">RUN</a>
    <a href="login.html">LOGIN</a>
    <a href="guide.html">GUIDE</a>
    <a href="entrez.html">ENTREZ</a>
    <a href="webeditor.php">EDIT</a>
    <a href="logout.php">LOGOUT</a>
    <button class="btn_change_theme" id="toggleBackgroundBtn">Change Theme</button>
</div>

<img src="figure/logo-transparent-png.png" width="360" class="imageContainer"/>
</br>
</br>
</br>
</br>
</br>
</br>
<div class="container">
    <div class="card">
    <h4>Server Status: <span id="status">...Connecting...</span></h4>
        <div class="clock" id="clock"></div>
        <hr>
        <h3>Upload PDB File</h3>
        <form action="upload.php" method="post" enctype="multipart/form-data">
            <div class="drop_box">
                <header>
                    <h4>Select Protein Strcture Here</h4>
                </header>
                <p>Files Supported: pdb</p>
                <input type="file" hidden accept=".pdb" id="fileID1" style="display:none;" name="file1" onchange="updateFileName('fileID1', 'fileName1')">
                <button class="btn" type="button" onclick="document.getElementById('fileID1').click()">Choose File</button>
                <span id="fileName1"></span>
            </div>
            </br>
            <h3>Protein Mutations</h3>
            <div class="drop_box">
                <header>
                    <h4>Record Your Mutations Here</h4>
                </header>
                <p>Seperate Category And Mutation With Tab(\t)</p>
                <textarea type="text" id="mutations" name="mutationstring" class="mutation-input" required placeholder="ASD     M1V"></textarea>
            </div>
            </br>
            <h3>Email</h3>
            <div class="drop_box">
                <header>
                    <h4>Give Us Your Email</h4>
                </header>
                <p>You'll Be Notified When Results Are Ready</p>
                <input type="text" id="searchStringAlphafold" name="searchString" class="search-input" required placeholder="spencerwwwwww@gmail.com">
            </div>
            </br>
            </br>
            <div class="drop_box_submit">
                <button class="btn_submit">Submit</button>
            </div>
        </form>
    </div>
    <br/>
    <br/>


</div>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<div class="footer">
    <!-- 底部内容 -->
    Department of Bioinformatics,
    Medical School of Soochow University <br>
    Contact: <a href="mailto:spencer-jrwang@foxmail.com">spencer-jrwang@foxmail.com</a><br>
    Source Code: <a href="https://github.com/Spencer-JRWang/APMA">Github</a>
</div>

<script>
    // 1. Show a welcome message when the homepage opens.
    window.onload = function() {
        alert('Welcome to SUDA APMA Server!');
        startClock();
    };

    // 2. Show a goodbye message when the homepage closes.
    window.onbeforeunload = function() {
        return 'Goodbye! See you again!';
    };

    // 3. Change the color of any paragraph text when it is clicked.
    document.querySelectorAll('p').forEach(function(p) {
        p.addEventListener('click', function() {
            p.style.color = p.style.color === 'blue' ? '#a3a3a3' : 'blue';
        });
    });

    // 4. Add a button to toggle the background theme of the homepage.
    document.getElementById('toggleBackgroundBtn').addEventListener('click', function() {
        document.body.classList.toggle('dark-theme');
    });

    // 5. Add a clock to the homepage that displays the current date and time.
    function startClock() {
        setInterval(() => {
            const now = new Date();
            const formattedTime = now.toLocaleString();
            document.getElementById('clock').textContent = formattedTime;
        });
    }

    function updateFileName(inputId, outputId) {
        var fileNameElement = document.getElementById(outputId);
        var fileInput = document.getElementById(inputId);
        var file = fileInput.files[0];
        fileNameElement.textContent = file ? file.name : '';
    }

    // 6. AJAX fetch status
    function fetchStatus() {
            fetch('upload.php')
                .then(response => response.json())
                .then(data => {
                    document.getElementById('status').innerText = data.status;
                })
                .catch(error => console.error('Error fetching status:', error));
        }

        setInterval(fetchStatus, 3000); // 每3秒更新一次状态
</script>
</body>
</html>


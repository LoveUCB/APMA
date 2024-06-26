<?php
// 定义文件路径
$countFile = 'figure/run_count.txt';
$lockFile = './run.lock';
// 设置目标文件路径
$uploadDir = '/home/wangjingran/APMA/data/';

// 创建计数器文件和队列文件（如果不存在）
if (!file_exists($countFile)) {
    file_put_contents($countFile, '0');
}


// 检查文件锁状态
$fp = fopen($lockFile, 'c+');
if (!flock($fp, LOCK_EX | LOCK_NB)) {
    $status = "...BUSY...";
} else {
    $status = "IDLE";
    // 释放文件锁
    flock($fp, LOCK_UN);
}
fclose($fp);

// 如果是上传请求
if ($_SERVER['REQUEST_METHOD'] === 'POST') {

    // 获取文件锁以防止并发运行
    $fp = fopen($lockFile, 'c+');

    // 如果无法获得文件锁，则继续排队
    if (!flock($fp, LOCK_EX | LOCK_NB)) {
        echo "...A process is already running. Please wait...";
        exit;
    }


    // 读取当前的计数器值
    $runCount = (int) file_get_contents($countFile);

    // 增加计数器值
    $runCount++;

    // 将新的计数器值写回文件
    file_put_contents($countFile, (string) $runCount);

    // 检查并处理上传文件错误
    function checkUploadError($file) {
        if ($file['error'] !== UPLOAD_ERR_OK) {
            echo "Error while uploading " . $file['name'];
            exit;
        }
    }

    // 移动上传文件到目标目录
    function moveUploadedFile($file, $uploadDir) {
        $fileName = basename($file['name']);
        $targetFile = $uploadDir . $fileName;
        if (move_uploaded_file($file['tmp_name'], $targetFile)) {
            echo "...{$fileName} uploading... success";
            echo "</br>";
        } else {
            echo "Error while uploading {$fileName}";
            echo "</br>";
            exit;
        }
    }

    // 检查文件上传错误
    checkUploadError($_FILES['file1']);

    

    // 处理上传的PDB文件
    moveUploadedFile($_FILES['file1'], $uploadDir);

   // Handle mutation information
    if (isset($_POST['mutationstring']) && !empty($_POST['mutationstring'])) {
        // Strip any \r (carriage return) characters from the input
        $mutation = str_replace("\r", "", $_POST['mutationstring']);
        $mutationFile = $uploadDir . 'position.txt';
        file_put_contents($mutationFile, $mutation . PHP_EOL, FILE_APPEND);
        echo "...recording mutations... success";
        echo "</br>";
    } else {
        echo "Invalid mutations";
        echo "</br>";
        exit;
    }

    // handle uniprot_id
    if (isset($_POST['uniprot_ID']) && !empty($_POST['uniprot_ID'])) {
        // Strip any \r (carriage return) characters from the input
        $uniprot_id = str_replace("\r", "", $_POST['uniprot_ID']);
        $mutationFile = $uploadDir . 'uniport_ID.txt';
        file_put_contents($mutationFile, $mutation . PHP_EOL, FILE_APPEND);
        echo "...recording UniProt ID... success";
        echo "</br>";
    } else {
        echo "Failed to record ID";
        echo "</br>";
        exit;
    }

    // 处理邮箱地址
    if (isset($_POST['searchString']) && filter_var($_POST['searchString'], FILTER_VALIDATE_EMAIL)) {
        $email = $_POST['searchString'];
        $emailFile = $uploadDir . 'email.txt';
        file_put_contents($emailFile, $email . PHP_EOL, FILE_APPEND);
        echo "...recording email... success";
        echo "</br>";
    } else {
        echo "Invalid email address";
        echo "</br>";
        exit;
    }

    // 运行Python脚本
    $output = shell_exec('nohup /home/wangjingran/miniconda3/envs/APMA/bin/python /home/wangjingran/APMA/Run.py > /home/wangjingran/APMA/run.log &');
    echo '<h2>...APMA is running now...</h2>';
    echo 'You will be notified when the process is finished';

} else {
    // 返回状态信息
    header('Content-Type: application/json');
    echo json_encode([
        'status' => $status
    ]);
}
?>

import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
import smtplib
from email.header import Header
from email.utils import formataddr
from email.mime.image import MIMEImage
from email.mime.multipart import MIMEMultipart
from email.message import EmailMessage
# from .. import __youremail__
# from .. import __emailkey__

def send_email(task_id, toEmailAddrs):
    # 设置服务器所需信息
    fromEmailAddr = 'spencer-jrwang@foxmail.com'
    password = '*****'
    #toEmailAddrs = ['3338561620@qq.com']
    
    # 设置email信息
    # ---------------------------发送带附件邮件-----------------------------
    # 邮件内容设置
    message =  MIMEMultipart()
    # 邮件主题
    message['Subject'] = 'Auto Protein Mutation Analyzer'
    # 发送方信息
    message['From'] = fromEmailAddr
    # 接受方信息
    message['To'] = toEmailAddrs[0]
    # 邮件正文内容
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Auto Protein Mutation Analyzer</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
        }}
        .container {{
            max-width: 600px;
            margin: 0 auto;
            padding: 20px;
            border: 1px solid #ccc;
        }}
        .signature {{
            margin-top: 20px;
            padding-top: 20px;
            border-top: 1px solid #ccc;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h2>Thank you for using deePheMut</h2>
        <p>Dear user,</p>
        <p>Your task has been successfully finished. Please visit:</p>
        <a href = "106.54.2.54/deePheMut/user_data/{task_id}" target="_blank">106.54.2.54/deePheMut/user_data/{task_id}</a>
        <p>Best,</p>
        
        <div class="signature">
            <p>Department of Bioinformatics</p>
            <p>Medical School of Soochow University</p>
            <p>You can send your feedback to <a href="mailto:spencer-jrwang@foxmail.com">spencer-jrwang@foxmail.com</a></p>
        </div>
    </div>
</body>
</html>
"""
    message.attach(MIMEText(html_content, 'html'))

    #message.attach(MIMEText('APMA analyzation is done.\n Please check the file attached bellow.\n--------------------------------------\nDeveloped by Spencer Wang', 'plain', 'utf-8'))
    
    # 构造附件
    # att_img2 = MIMEText(open(r'/home/wangjingran/APMA/Email/APMA_outcome.zip', 'rb').read(), 'base64', 'utf-8')
    # att_img2['Content-disposition'] = 'attachment;filename="APMA_outcome.zip"'
    # message.attach(att_img2)
    # ---------------------------------------------------------------------
    
    # 登录并发送邮件
    try:
        server = smtplib.SMTP('smtp.qq.com')  # qq邮箱服务器地址，端口默认为25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending outcome email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[INFO]error:", e)

def send_error_email(toEmailAddrs):
    # 设置服务器所需信息
    fromEmailAddr = 'spencer-jrwang@foxmail.com'  # 邮件发送方邮箱地址
    password = '*****'  # (注意不是邮箱密码，而是为授权码)
    message = EmailMessage()
    message['Subject'] = 'Auto Protein Mutation Analyzer'
    message['From'] = fromEmailAddr
    message['To'] = toEmailAddrs[0]
    message.set_content("Oops! Something went wrong, please check your files")
    with open('/home/wangjingran/APMA/run.txt', 'rb') as f:
        file_data = f.read()
        file_name = os.path.basename('/home/wangjingran/APMA/run.txt')
    message.add_attachment(file_data, maintype='text', subtype='plain', filename=file_name)


    try:
        server = smtplib.SMTP('smtp.qq.com')  # qq邮箱服务器地址，端口默认为25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending error email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[INFO]error:", e)

def send_start_email(task_id, toEmailAddrs):
    # 设置服务器所需信息
    fromEmailAddr = 'spencer-jrwang@foxmail.com'  # 邮件发送方邮箱地址
    password = '*****'  # (注意不是邮箱密码，而是为授权码)
    #toEmailAddrs = ['3338561620@qq.com']  # 邮件接受方邮箱地址，注意需要[]包裹，这意味着你可以写多个邮件地址群发
    
    # 设置email信息
    # ---------------------------发送带附件邮件-----------------------------
    # 邮件内容设置
    message =  MIMEMultipart()
    # 邮件主题
    message['Subject'] = 'Your Task has been successfully created'
    # 发送方信息
    message['From'] = fromEmailAddr
    # 接受方信息
    message['To'] = toEmailAddrs[0]
    # 邮件正文内容
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>deePheMut</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            line-height: 1.6;
        }}
        .container {{
            max-width: 600px;
            margin: 0 auto;
            padding: 20px;
            border: 1px solid #ccc;
        }}
        .signature {{
            margin-top: 20px;
            padding-top: 20px;
            border-top: 1px solid #ccc;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h2>Thank you for using deePheMut</h2>
        <p>Dear user,</p>
        <p>Your Task is running now, you can view the task status on:</p>
        <p> <a href = "106.54.2.54/deePheMut/user_data/{task_id}" target="_blank">106.54.2.54/deePheMut/user_data/{task_id}</a> </p>
        <p>Best,</p>
        
        <div class="signature">
            <p>Department of Bioinformatics</p>
            <p>Medical School of Soochow University</p>
            <p>You can send your feedback to <a href="mailto:spencer-jrwang@foxmail.com" target="_blank">spencer-jrwang@foxmail.com</a></p>
        </div>
    </div>
</body>
</html>
"""
    # 将HTML内容作为MIMEText的一部分添加到邮件中
    message.attach(MIMEText(html_content, 'html'))
 
    # 登录并发送邮件
    try:
        server = smtplib.SMTP('smtp.qq.com')  # qq邮箱服务器地址，端口默认为25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending start email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[ERROR] error:", e)




if __name__ == "__main__":
    send_start_email('219e7b7e-a071-4d5f-a424-190bce60878d',["3338561620@qq.com"])
    # send_email(["3338561620@qq.com"])
    # send_error_email(["3338561620@qq.com"])
    

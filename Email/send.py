import os
import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart
from email.header import Header
from email.utils import formataddr
from email.mime.image import MIMEImage
from email.message import EmailMessage

def send_email(task_id, toEmailAddrs):
    # Set server information
    fromEmailAddr = 'spencer-jrwang@foxmail.com'
    password = '*****'
    
    # Set email information
    # ---------------------------Send email with attachment-----------------------------
    # Email content settings
    message = MIMEMultipart()
    # Email subject
    message['Subject'] = 'Auto Protein Mutation Analyzer'
    # Sender information
    message['From'] = fromEmailAddr
    # Recipient information
    message['To'] = toEmailAddrs[0]
    # Email body content
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
        <h2>Thank you for using protPheMut</h2>
        <p>Dear user,</p>
        <p>Your task has been successfully finished. Please visit:</p>
        <a href="106.54.2.54/protPheMut/user_data/{task_id}" target="_blank">106.54.2.54/protPheMut/user_data/{task_id}</a>
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

    # Construct attachment
    # att_img2 = MIMEText(open(r'/home/wangjingran/protPheMut/Email/protPheMut_outcome.zip', 'rb').read(), 'base64', 'utf-8')
    # att_img2['Content-disposition'] = 'attachment;filename="protPheMut_outcome.zip"'
    # message.attach(att_img2)
    # ---------------------------------------------------------------------
    
    # Login and send email
    try:
        server = smtplib.SMTP('smtp.qq.com')  # QQ email server address, default port is 25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending outcome email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[INFO] error:", e)

def send_error_email(toEmailAddrs):
    # Set server information
    fromEmailAddr = 'spencer-jrwang@foxmail.com'  # Sender email address
    password = '*****'  # (Note: this is not the email password, but an authorization code)
    message = EmailMessage()
    message['Subject'] = 'protPheMut'
    message['From'] = fromEmailAddr
    message['To'] = toEmailAddrs[0]
    message.set_content("Oops! Something went wrong, please check your files")
    with open('/home/wangjingran/protPheMut/run.txt', 'rb') as f:
        file_data = f.read()
        file_name = os.path.basename('/home/wangjingran/protPheMut/run.txt')
    message.add_attachment(file_data, maintype='text', subtype='plain', filename=file_name)

    try:
        server = smtplib.SMTP('smtp.qq.com')  # QQ email server address, default port is 25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending error email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[INFO] error:", e)

def send_start_email(task_id, toEmailAddrs):
    # Set server information
    fromEmailAddr = 'spencer-jrwang@foxmail.com'  # Sender email address
    password = '*****'  # (Note: this is not the email password, but an authorization code)
    
    # Set email information
    # ---------------------------Send email with attachment-----------------------------
    # Email content settings
    message = MIMEMultipart()
    # Email subject
    message['Subject'] = 'Your Task has been successfully created'
    # Sender information
    message['From'] = fromEmailAddr
    # Recipient information
    message['To'] = toEmailAddrs[0]
    # Email body content
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>protPheMut</title>
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
        <h2>Thank you for using protPheMut</h2>
        <p>Dear user,</p>
        <p>Your Task is running now, you can view the task status on:</p>
        <p> <a href="106.54.2.54/protPheMut/user_data/{task_id}" target="_blank">106.54.2.54/protPheMut/user_data/{task_id}</a> </p>
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
    # Attach HTML content as part of the MIMEText
    message.attach(MIMEText(html_content, 'html'))
 
    # Login and send email
    try:
        server = smtplib.SMTP('smtp.qq.com')  # QQ email server address, default port is 25
        server.login(fromEmailAddr, password)
        server.sendmail(fromEmailAddr, toEmailAddrs, message.as_string())
        print('[INFO] sending start email success')
        server.quit()
    except smtplib.SMTPException as e:
        print("[ERROR] error:", e)

if __name__ == "__main__":
    send_start_email('219e7b7e-a071-4d5f-a424-190bce60878d', ["3338561620@qq.com"])
    # send_email(["3338561620@qq.com"])
    # send_error_email(["3338561620@qq.com"])

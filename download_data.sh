# Download data from a temporary FileSender archive containing the Mera MB dataset.
# The link is found by downloading it with a browser, then in Download, right-click on the data and copy download URL.
# See https://www.forbesconrad.com/blog/download-wetransfer-to-linux-server-wget/

wget --user-agent Mozilla/4.0 'https://filesender.renater.fr/download.php?token=6a412638-b939-4bd1-93eb-d91269af3ae5&files_ids=39705894' -O mera_data.tgz
tar xzvf mera_data.tgz

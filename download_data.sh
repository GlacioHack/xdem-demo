# Download data from a temporary FileSender archive containing the Mera MB dataset.
# The link is found by downloading it with a browser, then in Download, right-click on the data and copy download URL.
# See https://www.forbesconrad.com/blog/download-wetransfer-to-linux-server-wget/

wget  --user-agent Mozilla/4.0 'https://filesender.renater.fr/download.php?token=00e59d13-98dc-4b7f-9ffd-92bb28e06987&files_ids=43430563' -O mera_data.tgz
tar xzvf mera_data.tgz
rm mera_data.tgz

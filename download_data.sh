# Download data from a temporary FileSender archive containing the Mera MB dataset.
# The link is found by downloading it with a browser, then in Download, right-click on the data and copy download URL.
# See https://www.forbesconrad.com/blog/download-wetransfer-to-linux-server-wget/

wget --user-agent Mozilla/4.0 'https://filesender.renater.fr/download.php?token=bc73932e-4999-4237-8554-2e1a8862fcd2&files_ids=39956953' -O mera_data.tgz
tar xzvf mera_data.tgz
rm mera_data.tgz

# Download data from a temporary FileSender archive containing the Mera MB dataset.
# The link is found by downloading it with a browser, then in Download, right-click on the data and copy download URL.
# See https://www.forbesconrad.com/blog/download-wetransfer-to-linux-server-wget/

wget --user-agent Mozilla/4.0 'https://filesender.renater.fr/download.php?token=5aab3c13-cdfb-41f5-a930-8da6368e2659&archive_format=undefined&files_ids=18517548' -O mera_data.tgz
tar xzvf mera_data.tgz
mv xdem_mera_data/ psf_nepal_2022/data/mb_Mera
rm mera_data.tgz

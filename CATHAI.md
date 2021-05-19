# CATHAI 2.0
## Notes
* Need to write your own config.py! See the sample config.py included
* Uses shiny server to serve the R shiny apps that are proxied through the authenication dashboard portal
* Make sure that the shiny server is only bound on 127.0.0.1 or that the host has a firewall (either local system or NeCTAR security group e.g.) to block shiny server ports to outside sources

## Setup and running
Setup conda and python modules:
```
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
conda create -n cathai2 -c conda-forge python=3 uwsgi
conda activate cathai2
pip install -r requirements.txt

In R:
install.packages(c("flexdashboard", "visNetwork", "rjson", "DT", "shinyjs", "igraph"))
```

Run dev:
```
export FLASK_APP=app/__init__.py
flask run
```
or:
```
uwsgi --socket 0.0.0.0:5000 --protocol=http -w wsgi:app
```
or (only good for live reload testing; too slow otherwise; use uwsgi):
```
honcho start -e config.env -f Local
```

Run prod:
```
uwsgi --ini /var/www/cathai-dev/cathai2.ini
```

Setup app database and init admin user defined in config.py:
```
python manage.py recreate_db
python manage.py setup_prod
```

Setup nginx reverse proxy for python app, setup LE for prod:
```
apt update && apt install nginx
cp cathai-dev.beatsonlab.com.conf /etc/nginx/sites-available
ln -s /etc/nginx/sites-available/cathai-dev.beatsonlab.com.conf /etc/nginx/sites-enabled/cathai-dev.beatsonlab.com.conf
nginx -t && systemctl reload nginx
add-apt-repository ppa:certbot/certbot
apt install python-certbot-nginx
certbot -d cathai-dev.beatsonlab.com
```

Setup systemctl service for prod:
```
cp cathai2.service /etc/systemd/system/
systemctl enable cathai2.service
systemctl start cathai2.service
```

## Notes

* Redis and workers (via honcho) not currently setup. Using synchronous uwsgi

## Scope

* Search feature for metadata



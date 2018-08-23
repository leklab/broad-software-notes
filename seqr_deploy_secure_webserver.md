### Installing and deploying an https webserver for seqr
---
- By: Prech Brian Uapinyoying 
- Created: 07/26/2018
- Last updated: 08/22/2018
- <span style="color:red">**STILL INCOMPLETE**</span> - need to know how to secure the site with a real web cert

We have tested out Seqr on a test/development server, which is not robust or secure. To fix this problem, we will use gunicorn and nginx.

0. Background reading material:
	+ [Run a Django app with gunicorn in Ubuntu 16.04, Part I](http://rahmonov.me/posts/run-a-django-app-with-gunicorn-in-ubuntu-16-04/)
	+ [Run a Django app with gunicorn in Ubuntu 16.04, Part II](http://rahmonov.me/posts/run-a-django-app-with-nginx-and-gunicorn/)
	+ My seqr directory `${SEQR_INSTALL_DIR} = /home/<user>/seqr`

1. Install gunicorn - application server
```bash
sudo apt-get update && apt-get upgrade
sudo apt-get install gunicorn
```

2. Start a gunicorn http to test out, and close it out
```bash
# Make sure to shut down the Django test sever before-hand
gunicorn wsgi:application --keep-alive 5 --bind 10.4.34.3 # Substitute with the IP of your VM
# Ctr+C to shutdown
```

3. Create a self-signed certificate. This is NOT secure, but we will use for now.
```bash
cd ${SEQR_INSTALL_DIR}/code/seqr
openssl req -x509 -newkey rsa:2048 -keyout server.key -out server.crt -days 365 -nodes
```
	+ [What are self-signed certificates?](http://www.gerv.net/security/self-signed-certs/)
	+ [How to create a self signed certificate with openssl](https://stackoverflow.com/questions/10175812/how-to-create-a-self-signed-certificate-with-openssl)

4. Install nginx
```bash
sudo apt-get install nginx
```

5. Generate a socket file using gunicorn, then shutdown server again
```bash
gunicorn --bind unix:/home/<user>/seqr/code/seqr/seqr.sock wsgi:application
# Ctr+C to shutdown
```

6. Create the static files from the Seqr/Django app that Nginx needs
```bash
cd ${SEQR_INSTALL_DIR}/code/seqr
./manage.py collectstatic

# files generated in ${SEQR_INSTALL_DIR}/code/seqr/static
# '/home/<user>/seqr/code/seqr/static'
```
	+ [Thread about nginx static file serving confusion with root alias](https://stackoverflow.com/questions/10631933/nginx-static-file-serving-confusion-with-root-alias)

7. Generate a new nginx configuration file `/etc/nginx/sites-available/seqr`, will need sudo. Paste below inside 'seqr' config file and modify to fit your settings:
```bash
server {
	listen 443 ssl http2;   # See link below for info
    server_name 10.4.34.3;  # The IP

    # Specify where the self-signed certificates are that were generated from step 3
    ssl_certificate /home/<user>/seqr/code/seqr/server.crt;
    ssl_certificate_key /home/<user>/seqr/code/seqr/server.key;

    location = /favicon.ico { access_log off; log_not_found off; }

    # Specify where static files were generated from step 6
    location /static/ {
            root /home/<user>/seqr/code/seqr;
    }

    # Specify socket file generated from step 5
    location / {
            include proxy_params;
        proxy_pass http://unix:/home/<user>/seqr/code/seqr/seqr.sock;
    }
}
```
	+ [How to enable http2 for increased security on nginx](https://bjornjohansen.no/enable-http2-on-nginx)

8. Delete the default nginx configuration file from `/etc/nginx/sites-enabled/`
```bash
sudo rm /etc/nginx/sites-enabled/default
```

9. Create a soft-link of the seqr file you just made in `/etc/nginx/sites-available/` to `/etc/nginx/sites-enabled`
```bash
sudo ln -s /etc/nginx/sites-available/seqr /etc/nginx/sites-enabled/seqr
```

10. Reload the config file for nginx with the new settings
```bash
# sudo service nginx reload
sudo service nginx restart

# if the server was up from previous instances and wasnt shutdown
pgrep gunicorn
kill $(pgrep gunicorn)

# run again, can add --workers 3
gunicorn --bind unix:/home/<user>/seqr/code/seqr/seqr.sock wsgi:application
```

11. How to further [secure](https://bjornjohansen.no/securing-nginx-ssl) the site with a real certificate, then [optimize](https://bjornjohansen.no/optimizing-https-nginx) it.
    - To be continued...

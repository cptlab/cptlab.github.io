________________________   ______        ______       ___       __    ______  
__  ____/__  __ \__  __/   ___  / ______ ___  /_      __ |     / /_______  /_ 
_  /    __  /_/ /_  /      __  /  _  __ `/_  __ \     __ | /| / /_  _ \_  __ \
/ /___  _  ____/_  /       _  /___/ /_/ /_  /_/ /     __ |/ |/ / /  __/  /_/ /
\____/  /_/     /_/        /_____/\__,_/ /_.___/      ____/|__/  \___//_.___/ 
                                                                              

Model Deployment: Slade Matthews, July 2024
Site Design: Raymond Lui, 2023

The website is now fully cloned to github at https://github.com/cptlab/cptlab.github.io.

The root directory of the site is /srv on the server.

All flask apps are in /srv/apps/app[X] and call model.py in the local directory and specify an access port.

All model deployments are now listed on /srv/tools.html.

They can be linked to via <p><a href="/Model1">Model 1</a></p> made on the fly from /templates/index.html.

All html files in /templates determine the look of the page for individual tools.

The conda environment for running the apps is set in /etc/systemd/system/flask_app[X].service


       __ ___         
 |\ | /__  |  |\ | \/ 
 | \| \_| _|_ | \| /\ 
                      
The location and PORT (5000 etc) of the deployment is determined in:

/etc/nginx/conf.d/cptlab-web.conf (must match port specified in flask_app[X].py).

Here is where /ModelX is linked to the port where the app is running.
e.g. website/Model1 --> flask app running on port 5000

Remember to sudo systemctl restart nginx


  __          ___  _  _   _       
 /__ | | |\ |  |  /  / \ |_) |\ | 
 \_| |_| | \| _|_ \_ \_/ | \ | \| 
                                  
In order to make these deployments persistent (auto restart on boot) we have to run them via an /etc/services file which calls gunicorn (remember to change the port!)

Here is a list of tasks for making the deployment persistent:

0. Make a wsgi file in /apps/app[X] that calls the flask_app[X].py file (this is called in the .service file).
1. Edit cptlab-web.conf file as above
2. Create file in /etc/systemd/system such as /etc/systemd/system/flask_app1.service
3. Enable the new service sudo systemctl enable flask_app1.service
4. Reload the daemon sudo systemctl daemon-reload
4. Start the new service sudo systemctl restart flask_app1.service

The site is set up to run from /srv and this is set as "main" on cptlab.github.io
Once changes are made to the site you can send them to github via:
cd /srv/
git add .
git commit -m "Update readme"
git push origin main

Or if changes are made on Github you must:
git pull origin main

If there are changes on both sides they have to be resolved.

Github authentication is via personal access token, it must be renewed every so often.
cd /srv and do git remote -v to see the personal access tokens.

All apps are run from their apps directory and are linked to via tools.html

Note: I removed the DNS records pointing to Gituhub pages IP's 	185.199.108.153, (109, 110, 110) from the A record and removed the CNAME file from githup pages.
And removed the _github-challenge-cptlab-org txt file from the DNS.



├── apps
│   ├── app1
│   │   ├── flask_app1.py
│   │   ├── model.py
│   │   ├── __pycache__
│   │   │   ├── flask_app1.cpython-312.pyc
│   │   │   ├── model.cpython-311.pyc
│   │   │   ├── model.cpython-312.pyc
│   │   │   └── wsgi.cpython-312.pyc
│   │   ├── static
│   │   │   ├── cptlab_icon.png
│   │   │   └── logo.jpg
│   │   ├── templates
│   │   │   ├── error.html
│   │   │   ├── index.html
│   │   │   └── result.html
│   │   ├── test_model.py
│   │   └── wsgi.py
│   ├── app2
│   │   ├── flask_app2.py
│   │   ├── model.pt
│   │   ├── model.py
│   │   ├── __pycache__
│   │   │   ├── flask_app2.cpython-312.pyc
│   │   │   ├── model.cpython-311.pyc
│   │   │   ├── model.cpython-312.pyc
│   │   │   └── wsgi.cpython-312.pyc
│   │   ├── static
│   │   │   ├── cptlab_icon.png
│   │   │   └── logo.jpg
│   │   ├── templates
│   │   │   ├── index.html
│   │   │   └── result.html
│   │   ├── test_model.py
│   │   └── wsgi.py
│   ├── app3
│   │   ├── cleanup.py
│   │   ├── flask_app3.py
│   │   ├── model.py
│   │   ├── __pycache__
│   │   │   ├── flask_app3.cpython-312.pyc
│   │   │   ├── model.cpython-312.pyc
│   │   │   └── wsgi.cpython-312.pyc
│   │   ├── results
│   │   │   ├── Test_data_result.csv
│   │   │   └── Xu2012_result.csv
│   │   ├── static
│   │   │   ├── cptlab_icon.png
│   │   │   └── logo.jpg
│   │   ├── templates
│   │   │   ├── error.html
│   │   │   ├── index.html
│   │   │   └── result.html
│   │   ├── uploads
│   │   │   ├── Test_data.csv
│   │   │   └── Xu2012.csv
│   │   └── wsgi.py
│   ├── app4
│   │   ├── cleanup.py
│   │   ├── flask_app4.py
│   │   ├── model.py
│   │   ├── __pycache__
│   │   │   ├── flask_app2.cpython-312.pyc
│   │   │   ├── flask_app4.cpython-312.pyc
│   │   │   ├── model.cpython-311.pyc
│   │   │   ├── model.cpython-312.pyc
│   │   │   └── wsgi.cpython-312.pyc
│   │   ├── static
│   │   │   ├── cptlab_icon.png
│   │   │   ├── logo.jpg
│   │   │   └── temp
│   │   │       ├── 4-amino-2,6-dinitrotoluene.png
│   │   │       ├── aldosterone.png
│   │   │       ├── Benzaldehyde.png
│   │   │       ├──  benzenecarbaldehyde.png
│   │   │       ├── BMS-189453.png
│   │   │       ├── carbamate.png
│   │   │       ├── cyanocobalamin.png
│   │   │       ├── Lsd.png
│   │   │       ├── menthol.png
│   │   │       ├── nitrate.png
│   │   │       ├── N,N-Diethylaniline hydrochloride.png
│   │   │       ├── retinoic acid.png
│   │   │       ├── retinol.png
│   │   │       ├── semaglutide.png
│   │   │       ├── Spironolactone.png
│   │   │       └── sucrose.png
│   │   ├── templates
│   │   │   ├── error.html
│   │   │   ├── error.html.old
│   │   │   ├── index.html
│   │   │   └── result.html
│   │   ├── test_model.py
│   │   └── wsgi.py
│   ├── gunicorn_config.py
│   └── __pycache__
│       ├── flask_app1.cpython-311.pyc
│       └── gunicorn_config.cpython-311.pyc
├── CNAME
├── contact.html
├── css
│   └── styles.css
├── depREADME.txt
├── img
│   ├── cptlab_computer_aided_drug_design.png
│   ├── cptlab_high_performance_computing.png
│   ├── cptlab_icon.png
│   ├── cptlab_in_silico_toxicology.png
│   ├── cptlab_logo_black-transparent.png
│   ├── cptlab_logo.png
│   ├── cptlab_logo_white-transparent.png
│   ├── cptlab_machine_learning.png
│   ├── cptlab_molecular_modelling.png
│   ├── cptlab_slade_matthews.jfif
│   └── cptlab_translational_regulatory_science.png
├── index.html
├── js
│   ├── app.js
│   └── particles.js
├── publications.html
├── README.md
└── tools.html



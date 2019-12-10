N=$1
rm -rf _build _site
cd ../
jupyter-book build champsbook
cd champsbook
bundle exec jekyll serve --watch --port $N

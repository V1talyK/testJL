include .env
main: ver rt up ps push down

ci rt:
	time docker-compose build $@
push pull down ps:
	time docker-compose $@
run-ci:
	docker-compose run --rm -p 8888:8888 ci
run-rt:
	docker-compose run --rm -p 8888:8888 rt
log:
	docker-compose logs -f rt
up:
	docker-compose up -d rt
sh:
	docker exec -it `docker-compose ps -q` sh
ver:
	docker run --rm -it -v `pwd`/.env:/src/.env:cached $(REG)/$(GULP) ver

test:
	http :8888/test Cookie:ma_sql=30

status:
	docker-compose run --rm rt julia -e "Pkg.status()"

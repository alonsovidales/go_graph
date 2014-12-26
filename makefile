NO_COLOR=\033[0m
OK_COLOR=\033[32;01m
ERROR_COLOR=\033[31;01m
WARN_COLOR=\033[33;01m

format:
	@echo "$(OK_COLOR)==> Formatting$(NO_COLOR)"
	go fmt ./...

test: deps
	@echo "$(OK_COLOR)==> Testing$(NO_COLOR)"
	go test ./...

lint:
	@echo "$(OK_COLOR)==> Linting$(NO_COLOR)"
	golint .

all: format lint test

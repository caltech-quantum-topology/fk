import os
import yaml

if __name__ == "__main__":
    config_dir = "."

    configs = [] 

    for filename in os.listdir(config_dir):
        if filename.endswith(".yaml") or filename.endswith(".yml"):
            with open(os.path.join(config_dir, filename), 'r') as file:
                config = yaml.safe_load(file)
                configs.append(config)

    merged_config = {}
    for config in configs:
        if "name" not in config:
            continue
        if "braid" not in config:
            continue
        if "inversion" not in config:
            continue
        merged_config[config["name"]] = {}
        merged_config[config["name"]]["braid"] = config["braid"]
        merged_config[config["name"]]["inversion"] = config["inversion"]

    with open(os.path.join(config_dir, "all_knots.yaml"), 'w') as outfile:
        yaml.dump(merged_config, outfile, default_flow_style=None)

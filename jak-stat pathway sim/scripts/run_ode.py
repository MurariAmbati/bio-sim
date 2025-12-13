from jak_stat.model_ode import simulate_ode


def main():
    df, meta = simulate_ode(t_end=60.0, cytokine_level=1.0)
    print(meta)
    print(df.tail())


if __name__ == "__main__":
    main()

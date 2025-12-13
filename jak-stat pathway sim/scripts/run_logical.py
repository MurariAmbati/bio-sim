from jak_stat.model_logical import simulate_logical


def main():
    df, meta = simulate_logical(steps=30, cytokine_on=True, pulse=True, pulse_end_step=8)
    print(meta)
    print(df.tail())


if __name__ == "__main__":
    main()

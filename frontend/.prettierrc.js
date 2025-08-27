module.exports = {
    semi: true,
    trailingComma: 'all',
    singleQuote: true,
    printWidth: 100,
    tabWidth: 4,
    overrides: [
        {
            files: ['**/*.yml', '**/*.yaml'],
            options: {
                singleQuote: false,
                tabWidth: 2,
            },
        },
    ],
};

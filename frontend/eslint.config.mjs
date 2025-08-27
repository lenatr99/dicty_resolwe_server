import jsxA11y from 'eslint-plugin-jsx-a11y';
import tsEslint from 'typescript-eslint';
import importPlugin from 'eslint-plugin-import';
import reactPlugin from 'eslint-plugin-react';
import reactHooksPlugin from 'eslint-plugin-react-hooks';
import eslintPluginPrettierRecommended from 'eslint-plugin-prettier/recommended';
import globals from 'globals';
import eslintJs from '@eslint/js';
import { fixupPluginRules } from '@eslint/compat';

export default [
    {
        ignores: [
            '**/node_modules',
            '**/build',
            '**/config',
            '**/docker',
            '**/scripts',
            '**/.yarn',
            'package-lock.json',
            '.prettierrc.js',
            'public/config.js',
            'public/mockServiceWorker.js',
            'config',
            'scripts',
        ],
    },
    ...tsEslint.config(
        eslintJs.configs.recommended,
        ...tsEslint.configs.recommended,
        ...tsEslint.configs.strict,
        ...tsEslint.configs.recommendedTypeChecked,
    ),
    eslintPluginPrettierRecommended,
    jsxA11y.flatConfigs.recommended,
    {
        plugins: {
            import: fixupPluginRules(importPlugin),
            react: reactPlugin,
            'react-hooks': fixupPluginRules(reactHooksPlugin),
        },
        languageOptions: {
            globals: {
                ...globals.browser,
                ...globals.jest,
                require: true,
                process: true,
            },
            parser: tsEslint.parser,
            ecmaVersion: 2020,
            sourceType: 'module',
            parserOptions: {
                project: true,
                ecmaFeatures: {
                    jsx: true,
                },
            },
        },

        settings: {
            'import/resolver': {
                node: {
                    extensions: ['.js', '.jsx', '.ts', '.tsx', '.json'],
                    moduleDirectory: ['node_modules', 'src'],
                    paths: ['src'],
                },
            },
            react: {
                version: 'detect',
            },
        },

        rules: {
            'class-methods-use-this': 'off',
            'no-console': 'warn',
            'no-debugger': ['warn'],
            'no-nested-ternary': 'off',
            'no-param-reassign': [
                'error',
                {
                    props: false,
                },
            ],
            'no-restricted-globals': 'warn',
            'no-sequences': 'warn',
            'no-void': 'off',
            'operator-assignment': 'off',
            'prefer-destructuring': 'off',
            'spaced-comment': 'warn',

            'import/extensions': [
                'error',
                'never',
                {
                    svg: 'always',
                    woff: 'always',
                    woff2: 'always',
                    css: 'always',
                    png: 'always',
                    jpg: 'always',
                    ico: 'always',
                    json: 'always',
                },
            ],
            'import/no-duplicates': 'warn',
            'import/prefer-default-export': 'off',
            'import/order': 'warn',
            'import/newline-after-import': 'warn',
            'import/no-extraneous-dependencies': [
                'error',
                {
                    devDependencies: true,
                },
            ],

            'jsx-a11y/no-autofocus': 'off',
            'jsx-a11y/click-events-have-key-events': 'off',
            'jsx-a11y/no-static-element-interactions': 'off',

            'prettier/prettier': 'warn',

            'react/jsx-filename-extension': [
                'error',
                {
                    extensions: ['.tsx'],
                },
            ],
            'react/destructuring-assignment': 'off',
            'react/require-default-props': 'off',
            'react/no-unused-prop-types': 'warn',
            'react/no-danger': 'warn',
            'react/react-in-jsx-scope': 'off',
            'react/jsx-boolean-value': 'warn',
            'react/jsx-curly-brace-presence': 'warn',
            'react/jsx-props-no-spreading': 'warn',
            'react-hooks/exhaustive-deps': [
                'warn',
                {
                    additionalHooks:
                        '(useReport|useStateWithEffect|useUploadFileDrop|useFetchWithAbort|useFetchCallbackWithAbort|useUpdateEffect|useItemDrop)',
                },
            ],

            '@typescript-eslint/no-unsafe-member-access': 'off',
            '@typescript-eslint/no-unsafe-assignment': 'off',
            '@typescript-eslint/no-unsafe-call': 'off',
            '@typescript-eslint/no-unsafe-return': 'off',
            '@typescript-eslint/no-unsafe-argument': 'off',
            '@typescript-eslint/ban-ts-comment': 'off',
            '@typescript-eslint/restrict-template-expressions': [
                'error',
                {
                    allowAny: true,
                },
            ],
            '@typescript-eslint/explicit-module-boundary-types': 'off',
            '@typescript-eslint/no-dynamic-delete': 'off',
            '@typescript-eslint/no-unused-vars': [
                'warn',
                {
                    argsIgnorePattern: '^_',
                },
            ],
            '@typescript-eslint/no-unused-expressions': 'off',
            '@typescript-eslint/no-empty-function': 'off',
            '@typescript-eslint/no-empty-interface': 'warn',
            '@typescript-eslint/lines-between-class-members': 'off',
            '@typescript-eslint/no-misused-promises': [
                'error',
                {
                    checksVoidReturn: false,
                },
            ],
        },
    },
];

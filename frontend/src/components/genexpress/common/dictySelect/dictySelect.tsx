import { FormControl, InputLabel, Select, SelectProps } from '@mui/material';
import React, { ReactElement, ReactNode } from 'react';
import { generateRandomString } from 'utils/stringUtils';

type DictySelectProps = {
    children: ReactNode;
    label: string;
    value: unknown;
    disabled?: boolean;
    handleOnChange: SelectProps['onChange'];
    className?: string;
};

const DictySelect = ({
    children,
    label,
    value,
    handleOnChange,
    disabled,
    className,
}: DictySelectProps): ReactElement => {
    const labelId = `${generateRandomString(5)}Label`;

    return (
        <FormControl variant="outlined" className={className}>
            <InputLabel id={labelId}>{label}</InputLabel>
            <Select
                labelId={labelId}
                value={value}
                onChange={handleOnChange}
                label={label}
                disabled={disabled}
            >
                {children}
            </Select>
        </FormControl>
    );
};

export default DictySelect;

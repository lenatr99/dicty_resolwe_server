import React, { ReactNode, ReactElement, forwardRef } from 'react';
import { Tooltip } from '@mui/material';
import useForwardedRef from '../useForwardedRef';
import { StyledIconButton, StyledIconButtonProps } from './iconButtonWithTooltip.styles';

type IconButtonWithTooltipProps = {
    title: string;
    disabled?: boolean;
    children: ReactNode;
    onClick: (event: React.MouseEvent<HTMLButtonElement, MouseEvent>) => void;
} & Pick<StyledIconButtonProps, '$disablePadding'>;

const IconButtonWithTooltip = forwardRef<HTMLButtonElement, IconButtonWithTooltipProps>(
    (
        { title, disabled, children, $disablePadding, onClick }: IconButtonWithTooltipProps,
        ref,
    ): ReactElement => {
        const forwardedRef = useForwardedRef<HTMLButtonElement>(ref);

        return (
            <Tooltip title={title} aria-label="">
                <span>
                    <StyledIconButton
                        aria-label={title}
                        onClick={onClick}
                        disabled={disabled}
                        $disablePadding={$disablePadding}
                        ref={forwardedRef}
                    >
                        {children}
                    </StyledIconButton>
                </span>
            </Tooltip>
        );
    },
);

export default IconButtonWithTooltip;

import React, { ReactElement, useRef, useState } from 'react';
import { Menu, IconButton } from '@mui/material';
import { MoreVert as MoreIcon } from '@mui/icons-material';
import { MobileSection, DesktopSection, DictyToolbar } from './dictyAppBar.styles';

type AppBarProps = {
    desktopSection: ReactElement;
    mobileMenu?: ReactElement;
};

const DictyAppBar = ({ desktopSection, mobileMenu }: AppBarProps): ReactElement => {
    const mobileMoreButtonElement = useRef(null);
    const [isMobileMenuOpen, setIsMobileMenuOpen] = useState(false);

    const closeMobileMenu = (): void => {
        setIsMobileMenuOpen(false);
    };

    const openMobileMenu = (): void => {
        setIsMobileMenuOpen(true);
    };

    const renderMobileMenu = (
        <Menu
            anchorEl={mobileMoreButtonElement.current}
            anchorOrigin={{ vertical: 'top', horizontal: 'right' }}
            keepMounted
            transformOrigin={{ vertical: 'top', horizontal: 'right' }}
            open={isMobileMenuOpen}
            onClose={closeMobileMenu}
        >
            {mobileMenu}
        </Menu>
    );

    const alwaysDesktop = mobileMenu == null;

    return (
        <>
            <DictyToolbar disableGutters>
                <DesktopSection $alwaysVisible={alwaysDesktop}>{desktopSection}</DesktopSection>
                {!alwaysDesktop && (
                    <MobileSection>
                        <IconButton
                            ref={mobileMoreButtonElement}
                            aria-label="show more"
                            aria-haspopup="true"
                            onClick={openMobileMenu}
                            color="inherit"
                        >
                            <MoreIcon color="action" />
                        </IconButton>
                    </MobileSection>
                )}
            </DictyToolbar>
            {renderMobileMenu}
        </>
    );
};

export default DictyAppBar;

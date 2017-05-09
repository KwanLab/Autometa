package lu.uni.lcsb.vizbin;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.beans.PropertyVetoException;
import java.beans.VetoableChangeListener;
import java.beans.VetoableChangeSupport;

import org.apache.log4j.Logger;

public abstract class ObjectWithProperties {
	/**
	 * Default class logger.
	 */
	private static Logger								logger	= Logger.getLogger(ObjectWithProperties.class);

	/**
	 * Set of verification listeners called when some property is changed.
	 */
	private final VetoableChangeSupport	vchs		= new VetoableChangeSupport(this);

	/**
	 * Set of listeners called when some property is changed.
	 */
	private final PropertyChangeSupport	pchs		= new PropertyChangeSupport(this);

	public ObjectWithProperties() {
		// logger property listener
		PropertyChangeListener propertyChangeLogger = new PropertyChangeListener() {
			@Override
			public void propertyChange(final PropertyChangeEvent arg0) {
				logger.debug("Property changed: " + arg0.getPropertyName() + ". Old: " + arg0.getOldValue() + " New: " + arg0.getNewValue());
			}

		};
		addPropertyChangeListener(propertyChangeLogger);

	}

	/**
	 * Add property change listener that can VETO (cancel) the property change.
	 * 
	 * @param listener
	 *          property listener
	 */
	public final void addVetoablePropertyChangeListener(final VetoableChangeListener listener) {
		vchs.addVetoableChangeListener(listener);
	}

	/**
	 * Removes property change l;istener that can VETO from the available
	 * listeners.
	 * 
	 * @param listener
	 *          listeren to be removed
	 */
	public final void removeVetoablePropertyChangeListener(final VetoableChangeListener listener) {
		vchs.removeVetoableChangeListener(listener);
	}

	/**
	 * Adds property change listener that is fired after property is changed.
	 * 
	 * @param listener
	 *          listener to be added
	 */
	public final void addPropertyChangeListener(final PropertyChangeListener listener) {
		pchs.addPropertyChangeListener(listener);
	}

	/**
	 * Removes property change listener from the list of all available listeners.
	 * 
	 * @param listener
	 *          listener to be removed
	 */
	public final void removePropertyChangeListener(final PropertyChangeListener listener) {
		pchs.removePropertyChangeListener(listener);
	}

	/**
	 * Method that fires all vetoable property change listener for a given
	 * property.
	 * 
	 * @param propertyName
	 *          property that has changed
	 * @param oldValue
	 *          old value of the property
	 * @param newValue
	 *          new value of the property
	 * @throws PropertyVetoException
	 *           if the change shouldn't be made this exception is thrown
	 */
	protected final void fireVetoableChange(final String propertyName, final Object oldValue, final Object newValue) throws PropertyVetoException {
		vchs.fireVetoableChange(propertyName, oldValue, newValue);
	}

	/**
	 * Method fires property change listener for a given property.
	 * 
	 * @param propertyName
	 *          name of the property that has changed
	 * @param oldValue
	 *          old value of the property
	 * @param newValue
	 *          new value of the property
	 */
	protected final void firePropertyChange(final String propertyName, final Object oldValue, final Object newValue) {
		pchs.firePropertyChange(propertyName, oldValue, newValue);
	}

	/**
	 * Returns list of all property change listeners (includes only standard
	 * property change listeners, vetoable property change listeners are
	 * excluded).
	 * 
	 * @return list of property change listeners
	 */
	protected final PropertyChangeListener[] getPropertyChangeListeners() {
		return pchs.getPropertyChangeListeners();
	}

}
